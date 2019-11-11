import numpy as np

def DelinWatershed(queue, direction, d4 = [1,2,3,4], printflag = False):
	nx, ny = direction.shape
	
	#initilize a matrix to store the mask
	marked = np.zeros(direction.shape)
	
	#D4 neighbors
	kd = np.array([[0,-1],
					[-1,0],
					[0,1],
					[1,0]])
	'''
	#make masks of which cells drain down, up, left right
	down = np.zeros(direction.shape)
	up = np.zeros(direction.shape)
	left = np.zeros(direction.shape)
	right = np.zeros(direction.shape)
	
	down[direction==d4[0]]=1
	left[direction==d4[1]]=1
	up[direction==d4[2]]=1
	right[direction==d4[3]]=1
	'''
	
	nqueue = queue.shape[0]
	
	ii=1
	while nqueue>0 : 
		if printflag:
			print("lap "+str(ii)+" ncell "+str(nqueue))
		queue2=np.array([]).reshape(0,2)
		
		#loop through the queue
		for i in range(nqueue):
			xtemp=int(queue[i,0])
			ytemp=int(queue[i,1])
			#add one to the subbasin area for the summary
			marked[xtemp,ytemp]=1
			
			#look for cells that drain to this cell
			for d in range(4):
				xus=int(xtemp-kd[d,0])
				yus=int(ytemp-kd[d,1])
				if xus in range(nx) and yus in range(ny): #if in the domain bounds
					if not np.isnan(direction[xus,yus]) and marked[xus,yus]==0: #if in the mask
						if direction[xus,yus]==d4[d]: #if pointing to cell
							#print(c(xus,yus))
							marked[xus,yus]=1 #add the upstream cell to the mask
							queue2=np.vstack([queue2, np.array([[xus,yus]])]) # and then add the upstream cell to the queue
		if queue2.shape[0]>=1:
			queue=queue2.copy()
			nqueue=queue.shape[0]
			ii+=1
		else:
			nqueue=0
	'''
	xrange, yrange = np.where(marked==1)
	
	output_dict = {'watershed':marked,
					'xrange':xrange,
					'yrange':yrange}
	'''
	
	return marked
