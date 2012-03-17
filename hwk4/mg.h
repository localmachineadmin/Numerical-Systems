#define alpha 0.001
#define dt 0.005
#define noise 0.99
#define NPRE  3
#define NPOST 3
#define NGMAX 15

//method to determine which relaxer to use
int mode;
void addint(float ***uf, float ***uc, float ***res, int nf);
void copy(float ***aout, float ***ain, int n);
void fill0(float ***u, int n);
void interp(float ***uf, float ***uc, int nf);
void relax(float ***u, float ***rhs, int n);
void resid(float ***res, float ***u, float ***rhs, int n);
void rstrct(float ***uc, float ***uf, int nc);
void slvsml(float ***u, float ***rhs);
