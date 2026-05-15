#!/bin/bash

selfdir=$(dirname $0)

# =============================================================================
# Fault geometry — all zone lengths constant across mesh resolutions
# Grid center: x2 = (i2 - (N2-1)/2.0) * DX  (symmetric, x2=0 between cells)
# BC condition: abs(x2) > 40 km (velocity-strengthening, Dirichlet)
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

N2=32768 # n points on grid
DX=2.5  # cell size
Dc=0.001
WDIR=$selfdir/single_fault_benchmark/output_Dc001

if [ ! -e $WDIR ]; then

    echo adding directory "$WDIR"
    mkdir -p "$WDIR"
fi

# use --verbose=2 to output all parameters
# use --verbose=1 to output only parameters different from previous patch
OMP_NUM_THREADS=100 /home/alba/Documents/motorcycle/2d/antiplane/build/motorcycle-ap-ratestate-serial \
    --verbose 1 \
    --epsilon 1e-6 \
    --export-state \
    --export-stress \
    --export-netcdf \
    --export-netcdf-rate 20 \
    --export-netcdf-step 4 \
    --maximum-step 3.15e7 \
    --maximum-iterations 155000000 \
    --friction-law 1 <<EOF
# output directory
$WDIR
# Rigidity
30e3
# time interval
3e10
# number of faults
1
# grid dimension (N2)
$N2
# sampling (dx2)
$DX
#   n  tau0   mu0   sig   a   b   L   Vo   G/(2Vs)   Vl Dirichlet
$(echo "" | awk -v n2="$N2" -v dx="$DX" -v Dc="$Dc" '
    function abs(x){return (x>0)?x:-x};
    BEGIN{
        c=1;
        tau0_p=-1; mu0_p=-1; sig_p=-1;
        a_p=-1; b_p=-1; L_p=-1;
        Vo_p=-1; damping_p=-1;
        Vl_p=-1; dirichlet_p="T";
    }{
        for (i2=0; i2<n2; i2++) {
            x2=(i2-(n2-1)/2.0)*dx;
            tau0=-1;   # Initial shear traction (MPa); use -1 for steady-state
            L=Dc;      # Characteristic weakening distance (m)
            a=1e-2;    # Rate-dependent parameter (unitless)

            if (abs(x2) <= 2.5e3) {
                b = a + 4.0e-3; # Velocity-weakening zone (-2.5 km to +2.5 km)
            } else {
                b = a - 4.0e-3; # Velocity-strengthening zone and outer loading zones
            }

            mu0=0.6;   # Reference coefficient of friction (unitless)
            sig=1e2;   # Effective normal stress (MPa)
            Vo=1e-6;   # Reference slip-rate (m/s)
            Vl=1e-9;   # Loading rate (m/s)
            damping=5; # Radiation damping coefficient (MPa/m/s)

            if (abs(x2) > 40e3) {
                dirichlet = "T"; # Apply Dirichlet boundary condition
            } else {
                dirichlet = "F"; # Resolve friction law
            }

            # Check if all parameters are identical to the previous patch
            if ((tau0_p == tau0) && (mu0_p==mu0) && (sig_p==sig) && (a_p==a) && (b_p==b) &&
                (L_p==L) && (Vo_p==Vo) && (damping_p==damping) && (Vl_p==Vl) &&
                (dirichlet_p==dirichlet)) {
                # Print minus line number to save space (previous value is used)
                printf "%5d\n", -c;
            } else {
                # Print new set of parameters
                printf "%5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %d %10.2e %s\n",
                        c, tau0, mu0, sig, a, b, L, Vo, damping, Vl, dirichlet;
            }
            c++;
            tau0_p=tau0; mu0_p=mu0; sig_p=sig;
            a_p=a; b_p=b; L_p=L;
            Vo_p=Vo; damping_p=damping;
            Vl_p=Vl; dirichlet_p=dirichlet;
        }
    }')
# number of observation patches
3
# -----------------------------------------------------------------------------------
#   n fault     i2 rate
# -----------------------------------------------------------------------------------
    1     1 $(awk -v n2=$N2 'BEGIN{print int(n2/2)+1}') 1                     # Center of Fault 1
    2     1 $(awk -v n2=$N2 -v dx=$DX 'BEGIN{print int(n2/2)-(2400/dx)}') 1   # 2.4 km to the left of the center
    3     1 $(awk -v n2=$N2 -v dx=$DX 'BEGIN{print int(n2/2)+(2400/dx)}') 1   # 2.4 km to the right of the center
# number of events (not implemented)
0
EOF
