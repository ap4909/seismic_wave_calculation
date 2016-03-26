from scitools.std import *
import sys
infile=open('%s' %sys.argv[1],'r')
infile.readline()
lines=infile.readlines()
theta=float(raw_input('Receiver-model centre distance(degrees)='))		#user inputs position of receiver in model in degrees
omega=theta
theta=((theta)*2*pi)/360		#theta in model, converted to radians, as python works in radians
if theta> float(lines[0].split()[4])*2*pi/360 or theta<(-1*float(lines[0].split()[4])*2*pi/360):		#test to see if receiver is within model
		print 'Receiver distance outside model. Program now ending'
		sys.exit()
rp=float(raw_input('Ray Parameter(seconds/degree)='))			#user inputs ray parameter			
rp=rp/(2*pi/360)		#converting seconds per degree to seconds per radian
pl=raw_input('Plot?(Y/N)')
thetas=[]		
for line in lines:						#creates list of the epicentral distances(in radians) at which velocities are sampled in model from the edge of the model to zero degrees
	if float(line.split()[4])> 0:
		thetas.append((float(line.split()[4])*2*pi)/360)
	else:
		thetas.append((float(line.split()[4])*2*pi)/360)		
		break
length=len(thetas)		#used later for formation of depths list, and velocities list
for u in range(length-2,-1,-1):		#adds negative epicental distances to list
	thetas.append(-1*thetas[u])
depths=[]
depths.append(float(0))
for line in lines[length:]:				#list of depths in model of the layer boundaries. Ensures first depth is zero
	a=float(line.split()[2])	
	if a not in depths:
		depths.append(a)
radii=[]					
for line in lines:		#creates list of radii sampled by model
	if float(line.split()[3]) not in radii:
		radii.append(float(line.split()[3]))
b=range(len(lines)/length)			#number of sampled radii in model
velocities={}			
for index in b:
	velocities[index]={}			#dictionary with a number of sub-dictionaries the same as number of sample radii in model
c1=0		#layer counter
c2=0		#theta counter
for line in lines:					#velocities of each layer stored in separate sub dictionaries
	if float(line.split()[4])>0:
		velocities[c1][c2]=float(line.split()[7])
		c2+=1
	if float(line.split()[4])==0:				#centre of model reached
		velocities[c1][c2]=float(line.split()[7])
		c2+=1
		for u in range(length-2,-1,-1):		#adds the velocities already in sub dictionary to the sub dictionary, creating a symmetrical sampling of velocities, and covering the whole model from edge to edge			
			velocities[c1][c2]=velocities[c1][u] 	
			c2+=1	
		c1+=1		#counter layer moves down one
		c2=0		#theta counter set to zero; back at outer edge of model
infile.close()
z_vals=[]		#list to store depths of ray
z_vals.append(0)
deg=[]		#list to store epicentral distances of ray
deg.append(theta*360/(2*pi)) 
earth_rad=radii[0]	#radius of earth, as specified in model data
l=0		#layer counter
aver_r=(radii[0]+radii[1])/2	#average radius of a layer; used in calculation of i below
lower_r=radii[1]		#radii of the lower boundary of the layer that the ray is passing through; used in calculation of x below
r=radii[0]		#radius of ray; necessary for calculation of the turning depth
thickness=depths[1]-depths[0]	#initial thickness of layer
z=0		#initial depth
ttime=0		#travel time for the wave
eps=2**(-1*(24)) 	#machine epsilon; used for detection of turning point
layer=0		#used below; keeps track of whether the ray has remained in the same layer or not after most recent calculation of the path through a cell
step=1		#used to adjust layer counter(l) when ray passes into a new layer; also used in the calculation of the velocity through a cell; altered after turning point
for e in range(len(thetas)-1):			
		if theta<=thetas[e] and theta>thetas[e+1]:							#searches through thetas list to find where the receiver is and then sets the velocity to the corresponding cell 
			v=(velocities[l][e]+velocities[l][e+1]+velocities[l+1][e]+velocities[l+1][e+1])/4	#velocity is an average of the four corners of the cell
			break
i=asin((rp*v)/aver_r)			#first calculation of incident angle of ray
while True:				#infinite loop
	if rp*((velocities[l+1][e]+velocities[l+1][e+1])/2)/lower_r>=(1-eps) and step==1:		#test for the turning cell(the cell which contains the turning point) , using the values at the bottom of the cell			
		delta_r=(radii[l]-radii[l+1])				#thickness of cell
		vb=(velocities[l+1][e]+velocities[l+1][e+1])/2		#average velocity at the bottom of the cell, assuming a linear variation in the velocities within the cell
		vt=(velocities[l][e]+velocities[l][e+1])/2		#average velocity at the top of the cell
		rturn=((rp*(vb-vt)*radii[l]/delta_r)+rp*vt)/(1+(rp*(vb-vt)/delta_r))	#finding the turning radius assuming a linear variation in velocity with the cell
		if layer==1 or layer==2:		#If ray entered new layer as it entered the turning cell, velocity is set to the average velocity at the top of the cell
			v=vt
		if layer==0:		#if ray did not enter a new layer as it entered the turning cell, velocity is set according to the depth of the ray in the cell, assuming a linear relationship between velocity and depth
			v=vt+((vb-vt)/delta_r)*(z-depths[l])	
		i=asin((rp*v)/r)		#incident angle in turning cell
		y=(r-rturn)
		ttime+=(math.sqrt(y**2+(y*tan(i))**2)/v)*2		#travel time variable updated. Multiplied by two as at turning point ray as in the turning layer the ray is symmetric
		z+=(r-rturn)		#turning depth
		z_vals.append(z*-1)
		print z,'turning depth'
		d_theta=((r-rturn)*tan(i))/rturn
		theta_turn= theta-d_theta
		print 20-(theta_turn*360/(2*pi)),'turning theta'
		for i in range(2):		#ray is symmetric in the turning layer
			theta-=d_theta		
			deg.append(theta*360/(2*pi))
		z=depths[l]	
		z_vals.append(z*-1)
		step*=-1		#as ray is now heading upwards the layer counter will be now decrease by one each time a new layer is entered
		aver_r=(radii[l]+radii[l+step])/2
		thickness=(radii[l-1]-radii[l])/2
		for e in range(len(thetas)-1):		#searches through thetas list to find the correct e value			
			if theta<=thetas[e] and theta>thetas[e+1]:
				break							 	
	i=asin((rp*v)/aver_r)		#incident angle of ray
	x=(theta-thetas[e+1])*aver_r		#calculates width of a cell	
	idash = atan(x/thickness)	#"critical angle": the angle that determines if the ray passes into a cell to the right of it, or into a cell below it
	if i>idash:			#test to see if ray remains in same layer, and moves into an adjacent cell
		layer=0		
		ttime+=(math.sqrt(x**2 + (x/tan(i))**2)/v)		# travel time calculated assuming cell is a square and using pythagoras' theorem		
		theta=thetas[e+1]	
		e+=1		#ray is now in next cell
		if step==1:		#ray is going downward
			z+=x/tan(i)		#depth of ray entering next cell
			r-=x/tan(i)
		else:		#ray is going upward
			z-=x/tan(i)	
			r+=x/tan(i)		
		z_vals.append(z*-1)		#multiplied by -1 so that a plot gives a 'u' shaped curve
		deg.append((theta*360)/(2*pi))
		thickness=abs(depths[l+step]-z)			#new thickness of cell
		if e==len(thetas)-1:		#edge of model reached; infinite loop broken
			break
		v=(velocities[l][e]+velocities[l][e+1]+velocities[l+step][e]+velocities[l+step][e+1])/4		#next velocity calculated based on 4 corners of the cell
	if i<idash:		
		layer=1
		ttime+=(math.sqrt(thickness**2+(thickness*tan(i))**2)/v)
		f=thickness*tan(i)/(aver_r*(thetas[e+1]-thetas[e]))		#treating the cell as a square and finding the ratio of the horizontal component of the distance travelled by the ray to the average cell width
		d_theta=f*(thetas[e+1]-thetas[e])		#decrease in theta from when ray enters cell to where it enters the cell below it; using f compensates for the cell not being a square
		theta-=d_theta				
		deg.append((theta*360)/(2*pi))		
		l+=step				#altering the layer counter by one as moving into new layer
		z=depths[l]
		z_vals.append(z*-1)
		if l==0 and step==-1:		#surface of the earth has been reached, infinite loop is broken
			break
		r=radii[l]
		aver_r=(radii[l]+radii[l+step])/2
		if step==1:		#lower r only necessary when ray is going downward as used in calculation of turning point 		
			lower_r=radii[l+1]
		v=(velocities[l][e]+velocities[l][e+1]+velocities[l+step][e]+velocities[l+step][e+1])/4
		thickness= abs(depths[l+step]-depths[l])
	if i==idash:		#ray passes directly through a corner of the cell	
		layer=2
		l+=step
		e+=1
		theta=thetas[e]
		z=depths[l]
		z_vals.append(z*-1)
		thetas.append(theta)
		aver_r=(radii[l]+radii[l+step])/2
		if step==1:		
			lower_r=radii[l+1]
		r=radii[l]
		v=(velocities[l][e]+velocities[l][e+1]+velocities[l+step][e]+velocities[l+step][e+1])/4
		thickness=abs(depths[l+step]-depths[l])
if pl=='Y':
	plot(deg,z_vals)
print 20-(theta*360/(2*pi))
print 'Travel time(s)=%s' %ttime
print 'Theta(degrees)=%s'%str(theta*360/(2*pi))



	

