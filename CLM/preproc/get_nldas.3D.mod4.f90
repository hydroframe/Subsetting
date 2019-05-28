! Program to extract NLDAS forcing variables over a given problem domain
! and interpolate to problem grid via nearest neighbor

! Program reads user-specified batch file containing file paths and domain info
! Batch file is currently hard-wired as batch.get_nldas in current directory
! Batch file must be in the following format:
! "/path_to_top_level_nldas_dir/"
! "/path_to_output_directory/"
! "/path_to_lat-lon_file/lat-lon_file_name.txt"
! nx ny nz                      (note that for forcing data, nz always equals 1)
! x0 y0 z0                      (i.e., ComputationalGrid.Lower.X, --.Lower.Y, --.Lower.Z)
! dx dy dz                      (i.e., ComputationalGrid.DX, etc.)
! start_year start_month start_day start_hour (YYYY MM DD HH)
! end_year end_month end_day end_hour         (YYYY MM DD HH)

! Output vars:
! Downward SW radiation at surface (DSWR) [W/m^2]
! Downward LW radiation at surface (DLWR) [W/m^2]
! Precipitation (APCP) [kg/m^2 (accumulated per hour) --> convert --> kg/m^2/s (precipitation rate)]
! Air temperature (2m) (Temp) [K]
! U-wind speed (10m) (UGRD) [m/s]
! V-wind speed (10m) (VGRD) [m/s]
! Surface air pressure (Press) [Pa]
! Air specific humidity (2m) (SPFH) [kg/kg]

program extract_nldas

 ! Logicals for directory tests:
 logical :: dir_e

 ! Integers for array sizes, loops
 integer i,j,k,ni,nj,npdi,npdj,npdk,lwbdi,lwbdj,upbdi,upbdj,startyear,endyear,startmonth,endmonth
 integer startday,endday,t,t_day,t_month,t_year,t_hour,newmonth,starthour,endhour,done,interp
 integer nt,ntcount,tstart,tstop

 ! Reals for lat/lon, input values, output values, grid information
 real*8                dummy1,dummy2,x0,y0,x1,y1,z1,dx,dy,dz
 real*8                lat,lon,GLAT,GLON,GLAT1,GLON1,GLAT2,GLON2,a,b,c,d
 real*8,allocatable :: data_in(:,:,:),data_out(:,:,:),test(:,:)
 real*8,allocatable :: PD_lon(:,:),PD_lat(:,:),dat_in(:,:,:),dat_out(:,:,:)
 real*8,allocatable :: DSWR(:,:,:),DLWR(:,:,:),APCP(:,:,:),TMP(:,:,:),             &  
                       UGRD(:,:,:),VGRD(:,:,:),PRESS(:,:,:),SPFH(:,:,:),pmap(:,:,:)
 character*200         outdir,filename,gribfile,gribdate,gribparmextract
 character*200         timestart,timestop,timestep,filenumber,endfilenumber
 character*200         daydirectory,monthdirectory,yeardirectory,ntcheck
 character*200         domain_latlon,nldasdirectory

 ! Set done = 0
 done = 0

 ! Dataset constants:
 ! ni = number of longitude cells in NLDAS forcing dataset (constant)
 ! nj = number of latitude cells in NLDAS forcing dataset (constant)
 ! x0 = lowest longitude in NLDAS 8forcing dataset (constant)
 ! y0 = lowest latitude in NLDAS forcing dataset (constant)
 ni   =  464
 nj   =  224
 x0   = -124.938
 y0   =  25.063

 ! User input:
 print*, "READING BATCH FILE: batch.get_nldas"
 open(99,file="batch.get_nldas",status="unknown")
 read(99,*), nldasdirectory
 read(99,*), outdir
 read(99,*), domain_latlon
 read(99,*), npdi,npdj,nt
 read(99,*), x1,y1,z1
 read(99,*), dx,dy,dz
 read(99,*), startyear,startmonth,startday,starthour
 read(99,*), endyear,endmonth,endday,endhour
 read(99,*), interp
 close(99)

 ! Print start time, end time
 print*, "start date =",startyear,startmonth,startday,starthour
 print*, "end date =",endyear,endmonth,endday,endhour

 ! Create met forcing data if doesn't exist
 print*, "Check if output directory needed"
 inquire( file=trim(adjustl(outdir)), exist=dir_e )
 if (dir_e ) then
  print*, "Output directory exists"
 else
  print*, "Making output directory ("//trim(adjustl(outdir))//")"
  call system("mkdir "//trim(adjustl(outdir)))
 end if

 ! Set time to t_start
 print*, "Set up time slice and problem domain, allocate arrays"
 t_year=startyear
 t_month=startmonth
 t_day=startday
 t_hour=starthour
 t=1

 ! Allocate arrays
 allocate( PD_lon(npdi,npdj), PD_lat(npdi,npdj) )
 allocate( data_in(8,ni,nj), data_out(8,npdi,npdj) )
 allocate( DSWR(npdi,npdj,nt), DLWR(npdi,npdj,nt), APCP(npdi,npdj,nt),  TMP(npdi,npdj,nt), &
           UGRD(npdi,npdj,nt), VGRD(npdi,npdj,nt), PRESS(npdi,npdj,nt), SPFH(npdi,npdj,nt) )
 allocate( test(ni,nj) ) 

 ! Define formatting
 12345 Format (i6.6)
 98768 Format (i4.4)
 98767 Format (i4.4,i2.2)
 98766 Format (i4.4,i2.2,i2.2) 
 98765 Format (i4.4,i2.2,i2.2,i2.2)

 ! Read lat/lon values over problem domain
 print*, "Read in lat/lon over problem domain"
 open(20,file=trim(domain_latlon))
 do j=1,npdj
  do i=1,npdi
   ! read(20,*) PD_lon(i,j), PD_lat(i,j)
   read(20,*) PD_lat(i,j), PD_lon(i,j)
  enddo
 enddo
 close(20)

 ! Set lower/upper lat/lon bounds of problem domain
 ! (lower/upper i,j of NLDAS grid encompassing problem domain)
 print*, "Defining lower/upper lat/lon bounds of problem domain"
 lwbdi=int(((PD_lon(1,1)-x0)/0.125)-5)
 lwbdj=int(((PD_lat(1,1)-y0)/0.125)-5)
 upbdi=int(((PD_lon(npdi,npdj)-x0)/0.125)+5)
 upbdj=int(((PD_lat(npdi,npdj)-y0)/0.125)+5)
 lwbdi=5
 upbdj=450

 print*, "--> Latitude bounds of problem domain: ", PD_lat(1,1), PD_lat(npdi,npdj)
 print*, "--> Longitude bounds of problem domain:", PD_lon(1,1), PD_lon(npdi,npdj)
 print*, "--> Latitude index bounds of problem domain: ", lwbdj, upbdj
 print*, "--> Longitude index bounds of problem domain:", lwbdi, upbdi



 ! Get the NLDAS Cell map
 if (interp.eq.2) then
  allocate( pmap(npdi,npdj,6))
  print*, "Bilinear interpolation"
  do j=1,npdj                     !Loop over domain grid -- i,j
   do i=1,npdi
    lon   = PD_lon(i,j)              !Problem cell center
    lat   = PD_lat(i,j)              !Problem cell center  
    lwbdi=int(((lon-x0)/0.125)-5)
    lwbdj=int(((lat-y0)/0.125)-5)
    upbdi=int(((lon-x0)/0.125)+5)
    upbdj=int(((lat-y0)/0.125)+5)
       
    do jj=lwbdj,upbdj                  !Loop over NLDAS subgrid that overlies cell
     do ii=lwbdi,upbdi
      GLON1 = x0 + (dble(ii-1)*0.125)  !NLDAS grid cell center      (FOR FULL GRID LOOP)
      GLON2 = GLON1 + 0.125            !NEXT NLDAS grid cell center (one cell East of ii)
      GLAT1 = y0 + (dble(jj-1)*0.125)  !NLDAS grid cell center
      GLAT2 = GLAT1 + 0.125            !NEXT NLDAS grid cell center (one cell North of jj)

      if ( ((lat).ge.(GLAT1)).and.((lat).lt.(GLAT2)).and. &
           ((lon).ge.(GLON1)).and.((lon).lt.(GLON2)) ) then
         pmap(i,j,1)=GLON1
         pmap(i,j,2)=GLON2
         pmap(i,j,3)=GLAT1
         pmap(i,j,4)=GLAT2
         pmap(i,j,5)=ii
         pmap(i,j,6)=jj
            
      endif

       enddo ! ii
      enddo ! jj
     enddo ! i
    enddo ! j
   endif ! interp == 2



 ! Loop over timesteps:
 write(endfilenumber,98765) endyear,endmonth,endday,endhour
 do while (done.lt.1)

  print*, "******************************************************************************" 
  print*, "STARTING NT LOOP"
  print*, "******************************************************************************" 
  print*, " "

  ! Loop over nt files...
  tstart     = t
  do ntcount = 1,nt

   ! Marker
   print*, "******************************************************************************" 
   print*, "NTCOUNT:", ntcount
   print*, "DATE:   ", t_year,t_month,t_day,t_hour
   print*, "END:    ", endyear,endmonth,endday,endhour
   print*, " "

   ! Reset arrays to zero
   data_in = 0.d0; data_out = 0.d0

   ! Set directory and file names
   print*, "Set file name for GRIB input"
   write(filenumber,98765) t_year,t_month,t_day,t_hour
   write(daydirectory,98766) t_year,t_month,t_day
   write(monthdirectory,98767) t_year,t_month
   write(yeardirectory,98768) t_year

   ! Set grib file path, wgrib parameters
   ! If data in subdirectories by month:
   ! gribfile=trim(nldasdirectory)//"/"//trim(yeardirectory)//"/"//trim(monthdirectory)//"/" &
   !                              //trim(daydirectory)//"/"//trim(filenumber)//".nldasforce-a.grb"
   ! If data in subdirectories by year only:
   gribfile=trim(nldasdirectory)//"/"//trim(yeardirectory)//"/"//trim(daydirectory)//"/"//trim(filenumber)//".nldasforce-a.grb"
   gribparmextract=' | egrep "(:DSWRF:|:DLWRF:|:APCP:|:TMP:|:UGRD:|:VGRD:|:PRES:|:SPFH:)" | '

   ! Execute wgrib
   print*, "Execute wgrib (grb -> txt): "//trim(gribfile)//trim(gribparmextract)
   call system (" /data/LLNL/prov/wgrib/wgrib -s "//trim(adjustl(gribfile))//trim(adjustl(gribparmextract))//&
               " /data/LLNL/prov/wgrib/wgrib -i -s -text "//trim(gribfile)//" -o nldas_full.txt")

   ! Read data over full domain
   print*, "Read in data over full domain from nldas_full.txt"
   open(10,file="nldas_full.txt")
   do k=1,8                                   !Loop over vars
    read(10,*) dummy1, dummy2                 !(Skip first line -- nx,ny)
    do j=1,nj                                 !Loop over y
     do i=1,ni                                !Loop over x
      read(10,*) data_in(k,i,j)
     enddo
    enddo
   enddo
   close (10)

   ! MOVED OUTSIDE OF TIME LOOP -- VALUES ARE CONSTANT
   ! Read lat/lon values over problem domain
   ! print*, "Read in lat/lon over problem domain"
   ! open(20,file=trim(domain_latlon))
   ! do j=1,npdj                                !NOTE: Loops y-dimension first
   !  do i=1,npdi                               !      (must match code that created lat/lon file)
   !   ! read(20,*) PD_lon(i,j), PD_lat(i,j)    !NOTE: Reads lat,lon (NOT lon,lat)
   !     read(20,*) PD_lat(i,j), PD_lon(i,j)    !      (must also match code that created file)
   !  enddo
   ! enddo
   ! close(20)

   ! MOVED OUTSIDE OF TIME LOOP -- VALUES ARE CONSTANT 
   ! Set lower/upper lat/lon bounds of problem domain
   ! (lower/upper i,j of NLDAS grid encompassing problem domain)
   ! print*, "Define lower/upper lat/lon bounds of problem domain"
   ! lwbdi=int((PD_lon(1,1)-x0)/0.125-5)
   ! lwbdj=int((PD_lat(1,1)-y0)/0.125-5)
   ! upbdi=int((PD_lon(npdi,npdj)-x0)/0.125+5)
   ! upbdj=int((PD_lat(npdi,npdj)-y0)/0.125+5)
   ! print*, "-->Latitude index bounds of problem domain: ", lwbdj, upbdj
   ! print*, "-->Longitude index bounds of problem domain:", lwbdi, upbdi

   ! Select forcing data over problem domain
   ! Nearest neighbor algorithm
   if (interp.eq.1) then
    print*, "Retrieve met forcing for problem domain, interpolate via nearest neighbor"
    do j=1,npdj                          !Loop over domain grid -- i,j
     do i=1,npdi                    
      do jj=lwbdj,upbdj                  !Loop over NLDAS subgrid that overlies problem domain
       do ii=lwbdi,upbdi           
 
        GLAT = y0 + (dble(jj-1)*0.125)   !Latitude of NLDAS cell center   (FOR FULL GRID LOOP) 
        GLON = x0 + (dble(ii-1)*0.125)   !Longitude of NLDAS cell center   
        lat  = PD_lat(i,j)               !Latitude of domain cell center
        lon  = PD_lon(i,j)               !Longitude of domain cell center
 
        if ( ((lat).ge.(GLAT-0.0625)).and.((lat).lt.(GLAT+0.0625)).and.&
             ((lon).ge.(GLON-0.0625)).and.((lon).lt.(GLON+0.0625)) ) then      
 
           data_out(1:8,i,j) = data_in(1:8,ii,jj)          

        endif

       enddo 
      enddo
     enddo
    enddo
   endif !Endif @ nearest neighbor

   if (interp.eq.2) then
    print*, "Retrieve met forcing for problem domain, interpolate via bilinear interpolation"
    do j=1,npdj                     !Loop over domain grid -- i,j
     do i=1,npdi
      
         lon   = PD_lon(i,j)              !Problem cell center
         lat   = PD_lat(i,j)              !Problem cell center       
         GLON1=pmap(i,j,1)
         GLON2=pmap(i,j,2)
         GLAT1=pmap(i,j,3)
         GLAT2=pmap(i,j,4)
         ii=pmap(i,j,5)
         jj=pmap(i,j,6)

           do k=1,8
            a  = ((GLAT2-lat)/(GLAT2-GLAT1))*((GLON2-lon)/(GLON2-GLON1))*data_in(k,ii,jj)
            b  = ((GLAT2-lat)/(GLAT2-GLAT1))*((lon-GLON1)/(GLON2-GLON1))*data_in(k,ii+1,jj)
            c  = ((lat-GLAT1)/(GLAT2-GLAT1))*((GLON2-lon)/(GLON2-GLON1))*data_in(k,ii,jj+1)
            d  = ((lat-GLAT1)/(GLAT2-GLAT1))*((lon-GLON1)/(GLON2-GLON1))*data_in(k,ii+1,jj+1)
            data_out(k,i,j) = a+b+c+d
           enddo
     enddo ! i
    enddo ! j
   endif ! interp == 2

   !Separate variables from data_out into individual forcing vars:
   print*, "Separate forcing vars to individual arrays"
   do j=1,npdj
    do i=1,npdi
     TMP(i,j,ntcount)   = data_out(1,i,j)
     SPFH(i,j,ntcount)  = data_out(2,i,j)
     PRESS(i,j,ntcount) = data_out(3,i,j)
     UGRD(i,j,ntcount)  = data_out(4,i,j)
     VGRD(i,j,ntcount)  = data_out(5,i,j)
     DLWR(i,j,ntcount)  = data_out(6,i,j)
     APCP(i,j,ntcount)  = data_out(7,i,j) / 3600.                    ! conver kg/m^2/hr -> kg/m^2/s  (~ mm/s)
     DSWR(i,j,ntcount)  = data_out(8,i,j)
    enddo ! i
   enddo ! j

   ! Time Keeping
   ! Increment hour, day, year
   if (t_hour.eq.23) then

    ! Increment from 23 -> 0
    t_hour = 0

    ! See if need to increment month
    call months_end(t_year,t_month,t_day,newmonth)

    ! If not incrementing month, increment day only
    if (newmonth.eq.0) then
     t_day = t_day + 1

    ! If incrementing month, reset day -> 1, increment month
    else
     t_day = 1

     if (t_month.eq.12) then
      t_month = 1
      t_year  = t_year + 1
     else
      t_month = t_month + 1
     endif

    endif

   ! If t_hour != 23, increment hour only
   else
    t_hour = t_hour + 1
   endif
  
   ! Increment filenumber
   print*, "filenumber =", trim(filenumber),",  endfilenumber =",trim(endfilenumber)
   print*, " "
   t = t + 1

   ! Test if completed all timesteps
   if ((t_year.eq.endyear).and.(t_month.eq.endmonth).and.(t_day.eq.endday).and.(t_hour.eq.endhour)) then
    done = 1
   endif

  enddo ! ntcount
 
  ! Write forcing variables to PFB
  print*, "Write to pfb (pfb_write())"
  write(timestart,12345) tstart
  write(timestop,12345)  t-1

  filename=trim(adjustl(outdir))//"/NLDAS.Temp."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(TMP,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  filename=trim(adjustl(outdir))//"/NLDAS.SPFH."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(SPFH,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  filename=trim(adjustl(outdir))//"/NLDAS.Press."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(PRESS,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  filename=trim(adjustl(outdir))//"/NLDAS.UGRD."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(UGRD,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  filename=trim(adjustl(outdir))//"/NLDAS.VGRD."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(VGRD,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  filename=trim(adjustl(outdir))//"/NLDAS.DLWR."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(DLWR,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  filename=trim(adjustl(outdir))//"/NLDAS.APCP."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(APCP,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  filename=trim(adjustl(outdir))//"/NLDAS.DSWR."//trim(timestart)//"_to_"//trim(timestop)//".pfb"
  call pfb_write(DSWR,filename,npdi,npdj,nt,x1,y1,z1,dx,dy,dz)

  ! Reset arrays to zero
  DSWR = 0.d0; DLWR = 0.d0; APCP = 0.d0;
  TMP = 0.d0; UGRD = 0.d0; VGRD = 0.d0; PRESS = 0.d0

 enddo! (do while done.lt.1)
!print*, "Meteorological Forcing files are available in "//trim(adjustl(outdir))  

end program


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to handle days per month, leap years
subroutine months_end(year,mon,day,gotend)

 integer year,mon,day,gotend, month(12)

 ! February -- leap years
 if (MOD(year,4).eq.0) then
  month(2)=29
 else 
  if(MOD(year,400).eq.0) then
   month(2)=29
  else
   month(2)=28
  endif
 endif

 ! All other months
 month(1)=31
 month(3)=31
 month(4)=30
 month(5)=31
 month(6)=30
 month(7)=31
 month(8)=31
 month(9)=30
 month(10)=31
 month(11)=30
 month(12)=31

 ! Increment month
 if (day.eq.month(mon)) then
  gotend = 1 !Got to last day in the month
 else
  gotend=0
 endif 

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to write output as PFB
subroutine pfb_write(value,fname,nx,ny,nz,x1,y1,z1,dx,dy,dz)
  implicit none
  real*8        value(nx,ny,nz)
  real*8        dx,dy,dz,x1,y1,z1
  integer*4     i,j,k,nni,nnj,nnk,ix,iy,iz,                   &
                ns,rx,ry,rz,nx,ny,nz,nnx,nny,nnz,is
  character*200 fname

  ! ifort
  ! open(100,file=trim(fname),form='unformatted',accconvert='BIG_ENDIAN',status='unknown')

  ! gfortran
  open(100,file=trim(fname),form='unformatted',access='stream',convert='BIG_ENDIAN',status='unknown')

  nnx=nx; nny=ny; nnz=nz
  ns=1
  ix=0;iy=0;iz=0
  rx=0;ry=0;rz=0

  ! Start: writing of domain spatial information
  write(100) x1 !X
  write(100) y1 !Y 
  write(100) z1 !Z

  write(100) nx !NX
  write(100) ny !NY
  write(100) nz !NZ

  write(100) dx !DX
  write(100) dy !DY
  write(100) dz !DZ

  write(100) ns !num_subgrids
  ! End: writing of domain spatial information

  ! Start: loop over number of sub grids
  do is = 0, (ns-1)

   ! Start: writing of sub-grid spatial information
   write(100) ix
   write(100) iy
   write(100) iz
   write(100) nnx
   write(100) nny
   write(100) nnz
   write(100) rx
   write(100) ry
   write(100) rz
   ! End: writing of sub-grid spatial information

   ! Start: write in data from each individual subgrid
    !write(100) value((ix+1):(ix+nnx),(iy+1):(iy+nny),(iz+1):(iz+nnz))
    write(100) value
 ! do  k=iz +1 , iz + nnz
  !  do  j=iy +1 , iy + nny
  !   do  i=ix +1 , ix + nnx
  !    write(100) value(i,j,k)
  !   end do
  !  end do
  ! end do
   
   ! End: write data from each individual subgrid

  end do
  ! End: loop over number of sub grids

  close(100)
end subroutine



