(* xAct setup hook for emergent-metric/boundary calculations. *)
Needs["xAct`xTensor`"];
DefManifold[M,4,{a,b,c,d}];
DefMetric[-1,g[-a,-b],CD,PrintAs->"g"];
DefTensor[A[a],M];
DefTensor[th[],M];
