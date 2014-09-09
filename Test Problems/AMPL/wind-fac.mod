# wind-fac.mod   LOR2-AN-15-14
# created by Michal MICHNA, Polytechnika Gdanska,  June 2000.
#
# This file containes description of the winding factor of the electrical 
# machines.
# It is a test model for MINLP_BP solver - for Mixed Integer Nonlinear
# Programming problem
#

# sets 
# A set can be any unordered collection of objects pertinent to a model.

# parameters 
# A parameter is any numerical value pertinent to a model.

param  ms  := 3 integer;                # number of phases
param  p   := 2 integer;                # number of pole pairs
param  pi  := 3.141592654;              # pi
param  K   := 5 integer;                # harmonic order
param  ns  := 1;                        # coil span

# variabels
var q     >= 1, <=10, integer;  # number of slots per one phase and per one pole
var Nz    >= 1, integer;        # number of slots
var alfae := 1.5;
var tauz  := 1.0;               # slots pitch
var s     >= 1, integer;        # span  

# coil-group factor
var kz1;
var kz3;
var kz5;
var kz7;
var kz9;

# coil-span factor
var ks1;
var ks3;
var ks5;

# winding factor
var kw;
var kw1 >= 0.8;                 # winding factor for first harmonic
var kw3;
var kw5;


# objective function

minimize wind_fac:
        kw;

# constraints

subject to

   def_Nz:      Nz = 2 * ms * q * p;
   def_alfae:   alfae = (2 * pi * p)/Nz;
   def_tauz:    tauz = Nz / (2*p);
   def_s:       s = tauz - ns;
  
   def_kz1:     kz1 = sin(q * alfae / 2) / (q * sin( alfae / 2));
   def_ks1:     ks1 = sin((s * pi) / (tauz * 2));
   def_kw1:     kw1 = ks1 * kz1;

   def_kz3:     kz3 = sin(3 * q * alfae / 2) / (q * sin(3 * alfae / 2));
   def_ks3:     ks3 = sin((3 * s * pi) / (tauz * 2));
   def_kw3:     kw3 = ks3 * kz3;

   def_kz5:     kz5 = sin(5 * q * alfae / 2) / (q * sin(5 * alfae / 2));
   def_ks5:     ks5 = sin((5 * s * pi) / (tauz * 2));
   def_kw5:     kw5 = ks5 * kz5;


   def_kw:      kw = kw3 * kw3 + kw5 * kw5;

