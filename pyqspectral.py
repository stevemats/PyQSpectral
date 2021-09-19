#!/usr/bin/env/python
# https://github.com/stevemats/PyQSpectral

import sys,os,platform
import numpy as num
import matplotlib.pyplot as plt


def print_banner():
    """Check if python exist & print version"""

    #--------------------------------------------------------------------------
    # Check Py Version then Print main Banner
    #-------------------------------------------------------------------------- 
    try:
        impl = platform.python_implementation()
    except AttributeError:
        impl = "Python"

    version = platform.python_version()

    if '__pypy__' in sys.builtin_module_names:
        version += " (pypy %s)" % ".".join(str(v) for v in sys.pypy_version_info)

    try:
        which_python = os.path.relpath(sys.executable)
    except ValueError:

        """ On Windows, having  python executable is on a different drive
            than the sources cannot be relative.
        """
        which_python = sys.executable
    print('Python existence checker:', '=== %s %s (%s) ===' % (impl, version, which_python))
    sys.stdout.flush()
    
    print("""

    ██▓███ ▓██   ██▓  █████    ██████  ██▓███  ▓█████  ▄████▄  ▄▄▄█████▓ ██▀███   ▄▄▄       ██▓    
    ▓██░  ██▒▒██  ██▒▒██▓  ██▒▒██    ▒ ▓██░  ██▒▓█   ▀ ▒██▀ ▀█  ▓  ██▒ ▓▒▓██ ▒ ██▒▒████▄    ▓██▒    
    ▓██░ ██▓▒ ▒██ ██░▒██▒  ██░░ ▓██▄   ▓██░ ██▓▒▒███   ▒▓█    ▄ ▒ ▓██░ ▒░▓██ ░▄█ ▒▒██  ▀█▄  ▒██░    
    ▒██▄█▓▒ ▒ ░ ▐██▓░░██  █▀ ░  ▒   ██▒▒██▄█▓▒ ▒▒▓█  ▄ ▒▓▓▄ ▄██▒░ ▓██▓ ░ ▒██▀▀█▄  ░██▄▄▄▄██ ▒██░    
    ▒██▒ ░  ░ ░ ██▒▓░░▒███▒█▄ ▒██████▒▒▒██▒ ░  ░░▒████▒▒ ▓███▀ ░  ▒██▒ ░ ░██▓ ▒██▒ ▓█   ▓██▒░██████▒
    ▒▓▒░ ░  ░  ██▒▒▒ ░░ ▒▒░ ▒ ▒ ▒▓▒ ▒ ░▒▓▒░ ░  ░░░ ▒░ ░░ ░▒ ▒  ░  ▒ ░░   ░ ▒▓ ░▒▓░ ▒▒   ▓▒█░░ ▒░▓  ░
    ░▒ ░     ▓██ ░▒░  ░ ▒░  ░ ░ ░▒  ░ ░░▒ ░      ░ ░  ░  ░  ▒       ░      ░▒ ░ ▒░  ▒   ▒▒ ░░ ░ ▒  ░
    ░░       ▒ ▒ ░░     ░   ░ ░  ░  ░  ░░          ░   ░          ░        ░░   ░   ░   ▒     ░ ░   
            ░ ░         ░          ░              ░  ░░ ░                  ░           ░  ░    ░  ░
                                                                                
    [x] Modeler: Steve Matindi
    [x] Credit: Philip Mocz (2020), Princeton Univeristy
    [x] About: Schrodinger-Poisson system simulation using the spectral m3thod     
    [x] Usage Demo  : $ python pyqspectral.py                                                            
    """)

def main():
    #--------------------------------------------------------------------------
    # Simulation only on user Request
    #-------------------------------------------------------------------------- 
    menu_title = "Choose an option below to Continue:".title()
    print(menu_title)
    print('\n1. Run Simulation')
    print('2. Exit Simulation')
    while True:
        try:
            choice = int(input('Enter choice: '))
            if choice == 1:
                simulate()
                break
            elif choice == 2:
                print("\n", "GoodBye!")
                break
            else:
                print("Invalid choice. Enter a choice in menu. 1 or 2")
                main()
        except ValueError:
            print("Invalid choice. Enter 1 or 2")
    exit()

def simulate():
	""" Quantum simulation """
	
    #--------------------------------------------------------------------------
    # Simulation parameters
    #-------------------------------------------------------------------------- 
	N         = 512    # Spatial resolution
	t         = 0      # current time of the simulation
	tEnd      = 0.03    # time at which simulation ends
	dt        = 0.0001  # timestep
	tOut      = 0.0001  # draw frequency
	G         = 4000  # Gravitaitonal constant
	plotRealTime = True # switch on for plotting as the simulation goes along
	
    #--------------------------------------------------------------------------
    # Domain [0,1] x [0,1] \ Intial Condition
    #-------------------------------------------------------------------------- 
	L = 1    
	xlin = num.linspace(0,L, num=N+1)  # NB: x=0 & x=1 are the same point!
	xlin = xlin[0:N]                     # remove periodic point
	xx, yy = num.meshgrid(xlin, xlin)
	
	amp = 0.01
	sigma = 0.03
	rho = 0.9
	rho+= 2*amp*num.exp(-((xx-0.5)**2+(yy-0.5)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
	rho+= 1.5*amp*num.exp(-((xx-0.2)**2+(yy-0.7)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
	rho+= amp*num.exp(-((xx-0.4)**2+(yy-0.6)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
	rho+= amp*num.exp(-((xx-0.6)**2+(yy-0.8)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
	rho+= amp*num.exp(-((xx-0.8)**2+(yy-0.2)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
	rho+= amp*num.exp(-((xx-0.6)**2+(yy-0.7)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
	rho+= amp*num.exp(-((xx-0.7)**2+(yy-0.4)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
	rho+= amp*num.exp(-((xx-0.3)**2+(yy-0.3)**2)/2/sigma**2)/(sigma**3*num.sqrt(2*num.pi)**2)
    #--------------------------------------------------------------------------
    # normalize wavefunction to <|psi|^2>=1
    #-------------------------------------------------------------------------- 
	rhobar = num.mean( rho )
	rho /= rhobar
	psi = num.sqrt(rho)
	
	# Fourier Space Variables
	klin = 2.0 * num.pi / L * num.arange(-N/2, N/2)
	kx, ky = num.meshgrid(klin, klin)
	kx = num.fft.ifftshift(kx)
	ky = num.fft.ifftshift(ky)
	kSq = kx**2 + ky**2
	
	# Potential
	Vhat = -num.fft.fftn(4.0*num.pi*G*(num.abs(psi)**2-1.0)) / ( kSq  + (kSq==0))
	V = num.real(num.fft.ifftn(Vhat))
	
	# number of timesteps
	Nt = int(num.ceil(tEnd/dt))
	
	# prep figure
	fig = plt.figure(figsize=(6,4), dpi=80)
	grid = plt.GridSpec(1, 2, wspace=0.0, hspace=0.0)
	ax1 = plt.subplot(grid[0,0])
	ax2 = plt.subplot(grid[0,1])
	outputCount = 1
	
    #--------------------------------------------------------------------------
    # Simulation Main Loop
    #-------------------------------------------------------------------------- 
	for i in range(Nt):
		# (1/2) kick
		psi = num.exp(-1.j*dt/2.0*V) * psi
		
		# drift
		psihat = num.fft.fftn(psi)
		psihat = num.exp(dt * (-1.j*kSq/2.))  * psihat
		psi = num.fft.ifftn(psihat)
		
		# update potential
		Vhat = -num.fft.fftn(4.0*num.pi*G*(num.abs(psi)**2-1.0)) / ( kSq  + (kSq==0))
		V = num.real(num.fft.ifftn(Vhat))
		
		# (1/2) kick
		psi = num.exp(-1.j*dt/2.0*V) * psi
		
		# update time
		t += dt
		
		# plot in real time
		plotThisTurn = False
		if t + dt > outputCount*tOut:
			plotThisTurn = True
		if (plotRealTime and plotThisTurn) or (i == Nt-1):
			
			plt.sca(ax1)
			plt.cla()
			
			plt.imshow(num.log10(num.abs(psi)**2), cmap = 'inferno')
			plt.clim(-1, 2)
			ax1.get_xaxis().set_visible(False)
			ax1.get_yaxis().set_visible(False)	
			ax1.set_aspect('equal')	
			
			plt.sca(ax2)
			plt.cla()
			plt.imshow(num.angle(psi), cmap = 'bwr')
			plt.clim(-num.pi, num.pi)
			ax2.get_xaxis().set_visible(False)
			ax2.get_yaxis().set_visible(False)	
			ax2.set_aspect('equal')	
			
			plt.pause(0.001)
			outputCount += 1
			
			
    #--------------------------------------------------------------------------
    # Save figure \& Print banner on call
    #-------------------------------------------------------------------------- 
	plt.sca(ax1)
	plt.title(r'$\log_{10}(|\psi|^2)$')
	plt.sca(ax2)
	plt.title(r'${\rm angle}(\psi)$')
	plt.savefig('PyQspectral.png',dpi=240)
	plt.show()
	
	return 0
  

if __name__== "__main__":
  print_banner()
  main()
