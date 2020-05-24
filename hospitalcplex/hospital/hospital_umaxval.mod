int n = ...; // number of patients of a set N
int t = ...; // number of non-intersecting time slots of a set T
int m = ...; // number of operating rooms of a set M
int S = ...; // umaxvalue: amount of scenarios

int b_sj[s in 1..S][j in 1..n] = ...; // umaxvalue: time of the operation
int v_j[j in 1..n] = ...; // value of the operation
int r_j[j in 1..n] = ...; // ready time slot
int d_j[j in 1..n] = ...; // due time slot
int D_t[i in 1..t] = ...; // duration of time slot [minutes]
int delta_t[i in 1..t] = ...; // maximum time slot extension value [minutes]

{int} M_t[k in 1..t] = ...; // operating rooms which can be used simultaneously in time slots
{int} M_j[j in 1..n] = ...; // set of operating rooms in which operation j can be performed
{int} M_jt[j in 1..n][k in 1..t] = ...; // set of eligible operating rooms in time slot

{int} N_i[i in 1..m] = ...; // set of eligible operations for room i
{int} N_t[k in 1..t] = ...; // set of eligible operations for time slot
{int} N_it[i in 1..m][k in 1..t] = ...; // set of eligible operations for room i in time slot

{int} T_j[j in 1..n] = ...; // set of time slots eligible for operation j
{int} T_i[i in 1..m] = ...; // set of time slots in which room i is available
{int} T_ji[j in 1..n][i in 1..m] = ...; // set of time slots eligible for performing operation j in room i

float c_jit[j in 1..n][i in 1..m][k in 1..t] = ...; // cost
float e_it[i in 1..m][k in 1..t] = ...; // extension expense per time unit

dvar float z_it[i in 1..m][k in 1..t]; // the extension value of the time slot in operating room
dvar boolean x_jit[j in 1..n][i in 1..m][k in 1..t]; // patient is assigned to operating room in the time slot

// additional data (imported just for the results)
int situation = ...; // (n,t,m) combination index
int instance = ...; // instance index

float time;
execute {
  var before = new Date();
  time = before.getTime();
}

maximize 
  (sum(j in 1..n, k in T_j[j], i in M_jt[j][k])
    (v_j[j] - c_jit[j][i][k]) * x_jit[j][i][k]) -
  (sum(k in 1..t, i in M_t[k])
    e_it[i][k] * z_it[i][k]);
    
subject to {
   forall(s in 1..S, k in 1..t, i in M_t[k])
     sum(j in N_it[i][k]) b_sj[s][j] * x_jit[j][i][k] <= D_t[k] + z_it[i][k];
     
   forall(j in 1..n)
     sum(k in T_j[j], i in M_jt[j][k]) x_jit[j][i][k] <= 1;
     
   forall(k in 1..t, i in M_t[k])
     z_it[i][k] <= delta_t[k] * sum(j in N_it[i][k]) x_jit[j][i][k];
     
   forall(k in 1..t, i in M_t[k])
     sum(j in (asSet(1..n) diff N_it[i][k])) x_jit[j][i][k] == 0;
     
   forall(j in 1..n, k in 1..t, i in M_t[k])
     z_it[i][k] >= 0 && (x_jit[j][i][k] == 0 || x_jit[j][i][k] == 1);
}

execute {
  var after = new Date();
  time = after.getTime() - time;
}

int L_sit[s in 1..S][i in 1..m][k in 1..t] = sum(j in 1..n) x_jit[j][i][k] * b_sj[s][j];
int udl[s in 1..S] = sum(i in 1..m, k in 1..t) maxl(0, D_t[k] - L_sit[s][i][k]);

execute {
  var V_1 = cplex.getObjValue();
  var f = new IloOplOutputFile("..\\..\\hospitalcpp\\Data\\calcresults.inl", true);

  f.writeln("//(" + n + "," + t + "," + m + ") #" + instance);
  f.writeln("V_s[" + situation + "][" + instance + "] = " + V_1 + ";");
  f.writeln("time[" + situation + "][" + instance + "] = " + time + ";");
  for (var s = 0; s < S; s++) {
    f.writeln("udl[" + s + "][" + situation + "][" + instance + "] = " + udl[s+1] + ";");
  }

  f.writeln("");

  f.close();
}
