* FORM audit template: expand a symbolic two-loop cusp numerator.
* Run with: form cusp_trace.frm
V x1,y1,z1,w1,x2,y2,z2,w2;
Local n12 = w1*x1+w2*x2+y1*z1+y2*z2;
Local X12 = -(x1-x2)*(w1-w2)-(y1-y2)*(z1-z2);
Local I2 = -n12/(w1*x1*y1*z1*w2*x2*y2*z2*X12);
Print +s;
.end
