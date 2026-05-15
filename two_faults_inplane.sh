#!/bin/bash -e
# =============================================================================
# Fault geometry — all zone lengths constant across mesh resolutions
# Grid center: x2 = (i2 - (N2-1)/2.0) * DX  (symmetric, x2=0 between cells)
# BC condition: abs(x2) > 40 km (outer buffer, velocity-strengthening, Dirichlet)
#
#  Dc (m)  | DX (m) |   N2   | VW (m) | VS_L (m) | VS_R (m) | BC_L (m) | BC_R (m) | Total (m)
# ---------|--------|--------|--------|----------|----------|----------|----------|----------
#  0.0100  |   20   |   4096 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0080  |   20   |   4096 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0063  |   10   |   8192 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0050  |   10   |   8192 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0039  |   10   |   8192 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0030  |    5   |  16384 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0020  |    5   |  16384 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0014  |  2.5   |  32768 |  5000  |  37500   |  37500   |   960    |   960    |  81920
#  0.0010  |  2.5   |  32768 |  5000  |  37500   |  37500   |   960    |   960    |  81920
# =============================================================================
selfdir=$(dirname $0)

N2=16384 # n points on grid
DX=5  # cell size (m)
Dc=0.002

# fault spacing: 0.2e3 0.5e3 1e3 1.7e3 3e3 5e3 10e3 20e3 30e3 60e3 120e3

for dist in 0.2e3 0.5e3 1e3 1.7e3 3e3 5e3 10e3 20e3 30e3 60e3 120e3; do
    WDIR=$selfdir/inplane_two_fault/002/output_Dc002_${dist%e3}km

	if [ ! -e $WDIR ]; then
		echo adding directory "$WDIR"
		mkdir -p "$WDIR"
	fi
	
	# use --verbose=2 to output all parameters
	# use --verbose=1 to output only parameters different from previous patch
	OMP_NUM_THREADS=100 /home/alba/Documents/motorcycle/2d/planestrain/build/motorcycle-ps-ratestate-serial \
		--verbose 1 \
		--verbose 1 \
		--epsilon 1e-6 \
		--export-state \
		--export-stress \
		--export-netcdf \
		--export-netcdf-rate 20 \
		--export-netcdf-step 4 \
		--maximum-step 3.15e7 \
		--maximum-iterations 100500000 \
		--friction-law 1 <<EOF
# output directory
$WDIR
# Lame parameters (lambda, mu)
30e3 30e3
# time interval
3.5e10
# number of faults
2
# grid dimension (N2)
$N2
# sampling (dx2)
$DX
#   n  tau0   mu0   sig   a   b   L   Vo   G/(2Vs)   Vl Dirichlet
$(echo "" | awk -v n2="$N2" -v dx="$DX" -v Dc="$Dc" '
	function abs(x){return (x>0)?x:-x};
	function max(x,y){return (x>y)?x:y};
	function boxcar(x){return (x>=-0.5 && x<=0.5)?1:0};
	function heavi(x){return (x>0)?1:0};
	function ramp(x){return x*boxcar(x-0.5)+heavi(x-1)};
	function asinh(x){return log(x+sqrt(1+x^2))};
	function sinh(x){return (exp(x)-exp(-x))/2};
	BEGIN{
	c=1;
	# initial parameters for reuse
	tau0_p=-1;mu0_p =-1;sig_p=-1;
	a_p=-1;b_p=-1;L_p=-1;
	Vo_p=-1;damping_p=-1;
	Vl_p=-1;dirichlet_p="T";
	}{
	    for (i2=0; i2<n2; i2++) {
            x2=(i2-(n2-1)/2.0)*dx; # symmetric grid: x2=0 between the two central cells
            tau0=-1;   # Initial shear traction (MPa); use -1 for steady-state
            L=Dc;      # Characteristic weakening distance (m)
            a=1e-2;    # Rate-dependent parameter (unitless)

            if (abs(x2) <= 2.5e3) {
                b = a + 4.0e-3; # Velocity-weakening zone (5 km)
            } else if (abs(x2) <= 40e3) {
                b = a - 4.0e-3; # Velocity-strengthening zone (37.5 km each side)
            } else {
                b = a - 4.0e-3; # Outer buffer (960 m each side, Dirichlet BC)
            }

		# reference coefficient of friction (unitless)
		mu0=0.6;
		# effective normal stress (MPa)
		sig=1e2;
		# reference slip-rate (m/s)
		Vo=1e-6;
		# loading rate (m/s)
		Vl=1e-9;
		# radiation damping coefficient ( G/(2 Vs) in units of MPa/m*s )
		damping=5;
		if (abs(x2) > 40e3){
			# apply Dirichlet boundary condition over entire outer buffer
			dirichlet="T"
		} else {
			# resolve rate and state dependence of friction
			dirichlet="F";
		}
		if (1==boxcar((x2-0.5e3)/1e3)){
			# initial stress triggers nucleation (MPa)
			tau0=(mu0+(a-b)*log(1e-10/Vo))*sig;
		}
		# check if all parameters are identical to that of the previous patch
		if ((tau0_p == tau0) && (mu0_p==mu0) && (sig_p==sig) && (a_p==a) && (b_p==b) && (L_p==L) && 
		    (Vo_p==Vo) && (damping_p==damping) && (Vl_p==Vl) && (dirichlet_p==dirichlet)){
			# print minus line number to save space (previous value is used)
			printf "%5d\n",-c;
		} else {
			# print new set of parameters
			printf "%5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %d %10.2e %s\n", 
	        			c,  tau0,   mu0,   sig,     a,     b,     L,    Vo, damping, Vl, dirichlet;
		}
		c++;
		# save current parameters for reuse
		tau0_p=tau0;mu0_p=mu0;sig_p=sig;
		a_p=a;b_p=b;L_p=L;
		Vo_p=Vo;damping_p=damping;
		Vl_p=Vl;dirichlet_p=dirichlet;
	}
}')
# distance of second fault
$dist
#   n  tau0   mu0   sig   a   b   L   Vo   G/(2Vs)   Vl Dirichlet
$(echo "" | awk -v n2="$N2" -v dx="$DX" -v Dc=$Dc '
	function abs(x){return (x>0)?x:-x};
	function max(x,y){return (x>y)?x:y};
	function boxcar(x){return (x>=-0.5 && x<=0.5)?1:0};
	function heavi(x){return (x>0)?1:0};
	function ramp(x){return x*boxcar(x-0.5)+heavi(x-1)};
	function asinh(x){return log(x+sqrt(1+x^2))};
	function sinh(x){return (exp(x)-exp(-x))/2};
	BEGIN{
	c=1;
	# initial parameters for reuse
	tau0_p=-1;mu0_p =-1;sig_p=-1;
	a_p=-1;b_p=-1;L_p=-1;
	Vo_p=-1;damping_p=-1;
	Vl_p=-1;dirichlet_p="T";
	}{
	    for (i2=0; i2<n2; i2++) {
            x2=(i2-(n2-1)/2.0)*dx; # symmetric grid: x2=0 between the two central cells
            tau0=-1;   # Initial shear traction (MPa); use -1 for steady-state
            L=Dc;      # Characteristic weakening distance (m)
            a=1e-2;

            if (abs(x2) <= 2.5e3) {
                b = a + 4.0e-3; # Velocity-weakening zone (5 km)
            } else if (abs(x2) <= 40e3) {
                b = a - 4.0e-3; # Velocity-strengthening zone (37.5 km each side)
            } else {
                b = a - 4.0e-3; # Outer buffer (960 m each side, Dirichlet BC)
            }

		# reference coefficient of friction (unitless)
		mu0=0.6;
		# effective normal stress (MPa)
		sig=1e2;
		# reference slip-rate (m/s)
		Vo=1e-6;
		# loading rate (m/s)
		Vl=1e-9;
		# radiation damping coefficient ( G/(2 Vs) in units of MPa/m*s )
		damping=5;
		if (abs(x2) > 40e3){
			# apply Dirichlet boundary condition over entire outer buffer
			dirichlet="T"
		} else {
			# resolve rate and state dependence of friction
			dirichlet="F";
		}
		if ((tau0_p == tau0) && (mu0_p==mu0) && (sig_p==sig) && (a_p==a) && (b_p==b) && (L_p==L) && 
		    (Vo_p==Vo) && (damping_p==damping) && (Vl_p==Vl) && (dirichlet_p==dirichlet)){
			# print minus line number to save space (previous value is used)
			printf "%5d\n",-c;
		} else {
			# print new set of parameters
			printf "%5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %d %10.2e %s\n", 
	        			c,  tau0,   mu0,   sig,     a,     b,     L,    Vo, damping, Vl, dirichlet;
		}
		c++;
		tau0_p=tau0;mu0_p=mu0;sig_p=sig;
		a_p=a;b_p=b;L_p=L;
		Vo_p=Vo;damping_p=damping;
		Vl_p=Vl;dirichlet_p=dirichlet;
	}
}')
# number of observation patches
6
# -----------------------------------------------------------------------------------
#   n fault     i2 rate
# -----------------------------------------------------------------------------------
    1     1 $(awk -v n2=$N2 'BEGIN{print int(n2/2)+1}') 1                     # Center of Fault 1
    2     2 $(awk -v n2=$N2 'BEGIN{print int(n2/2)+1}') 1                     # Center of Fault 2
    3     1 $(awk -v n2=$N2 -v dx=$DX 'BEGIN{print int(n2/2)-(2400/dx)}') 1   # 2.3 km to the left of the center for Fault 1
    4     1 $(awk -v n2=$N2 -v dx=$DX 'BEGIN{print int(n2/2)+(2400/dx)}') 1   # 2.3 km to the right of the center for Fault 1
    5     2 $(awk -v n2=$N2 -v dx=$DX 'BEGIN{print int(n2/2)-(2400/dx)}') 1   # 2.3 km to the left of the center for Fault 2
    6     2 $(awk -v n2=$N2 -v dx=$DX 'BEGIN{print int(n2/2)+(2400/dx)}') 1   # 2.3 km to the right of the center for Fault 2
# number of events (not implemented)
0
EOF
done