from geocusp.continuum import run_entropy_transport
from geocusp.export import export_1d_field, export_openfoam_table
r=run_entropy_transport()
field={'sigma_IR':r.sigma,'m_eff':0.7*(1+0.8*__import__('numpy').exp(-(r.x/(0.7*4.0))**2))}
print(export_1d_field(r.x,field)); print(export_openfoam_table(r.x,r.sigma))
