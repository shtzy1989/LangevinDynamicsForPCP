#ifndef MD_CONST_H
#define MD_CONST_H

const double SIGMA2RMIN     = 1.12246204830937; // 2^(1/6)
const double PI             = 3.1415926535897932384626433832795029;
const double PI2            = PI * 2.0;
const double PIHALF         = 1.5707963267949;
const double PIR            = 0.318309886183791;
const double PIROOT         = 1.77245385090552;
const double PIROOT2        = 3.54490770181103;
const double PIROOTR        = 0.564189583547756;
const double RAD2DEG        = 57.29577951308232088; // 180 / PI;
const double DEG2RAD        = 1.0 / RAD2DEG; // PI / 180;
const double P1             = 0.0730f;
const double P2             = 0.9212;
const double P3             = 6.2105;
const double P4             = 15.2258;
const double P5             = 1.2535;
const double d_offset       = -0.09;
const double VELOCITY       = 20.45482706;
const double AMBERCHARGE    = 18.2223;
const double MINCUT         = 0.00001;
const double KB             = 1.38e-23;
const double PLANCK         = 6.626e-34;
const double RPLANCK        = 1.054571800e-34;
const double NA             = 6.02214199e+23;
const double GASCONST       = 8.314;
const double E2COULOMB      = 1.602188e-19;
const double C2STATC        = 3.00e9;
const double CHARGECONSTANT = 332.063709667747;
const double SMALL          = 1.0e-9;
const double HALF           = 0.5;
const double THIRD          = 0.333333333333333333;
const double QUATER         = 0.25;

const double MPS2APFS       = 1.0e-5;
const double APFS2MPS       = 1.0e5;

const double MPS22APFS2     = 1.0e-20;
const double APFS22MPS2     = 1.0e20;

const double BOHR2ANGSTROM  = 0.5291772109038;
const double HARTREE2J      = 4.3597447222071e-18;
const double ELECTRONVOLT2J = 1.602176634e-19;
const double RYDBERG2J      = HARTREE2J / 2.0;
const double RYDBERG2EV     = 13.605693122994;

#endif
