# script that writes co_chemsource subroutine from yaml file
# no custom orders are implemented
# by default considers reactions as reversible, is forward you need backward rate to be 0
# created with no love and absolutely no testing
# use at your own peril
# m.fabiani, 05/05/2025 ei fu siccome immobile

# adapt to work into FLINT workspace
# add CLI argument for mechanism name
# refactor to properly manage Lindemann and Troe falloff reactions
# tested on the already present mechanisms
# m.grossi, summer 2025

import cantera as ct
import sys

def find_indices(lst, condition):
   return [i for i, elem in enumerate(lst) if condition(elem)]

def stoich_to_mult(var: str, n: int) -> str:
    """Convert exponentiation to repeated multiplication or return '1.0' if zero power."""
    if n <= 0:
        return '1.0'
    elif n == 1:
        return var
    else:
        return '*'.join([var] * int(n))
    
def format_eff(val: float) -> str:
    """Simplify coefficient formatting."""
    if abs(val - 1.0) < 1e-8:
        return ''  # omit '*1.0'
    else:
        return f'*{val:.8g}'  # truncate trailing zeros
    
# Main script

print()
print( 'FLINT-Kinetix ' )
print( 'created with no love and partial testing')
print()

# Input file definition
if len(sys.argv) < 2:
  print("Usage: python YTF.py <mech>")
  sys.exit(1)

mech = sys.argv[1]

infilename = mech+".yaml"
outfilename = mech+".f90"
mechname = mech

all_species = ct.Species.list_from_file(infilename)
gas = ct.Solution(infilename)
ref_phase = ct.Solution(thermo='ideal-gas', kinetics='gas', species=all_species)
all_reactions = ct.Reaction.list_from_file(infilename, ref_phase)
species = gas.species_names
nr = len(all_reactions)

# Start writing the Fortran file

f = open(outfilename, "w")
f.write("module "+mechname+"_mod\n")
f.write("implicit none\n")
f.write("contains\n")
f.write("subroutine "+mechname+"(roi,temp,omegadot)\n")
f.write("use FLINT_Lib_Thermodynamic\n")
f.write("use FLINT_Lib_Chemistry_data\n")
f.write("use FLINT_Lib_Chemistry_falloff\n")
f.write("implicit none\n")
f.write("real(8), intent(inout)  :: roi(ns)\n")
f.write("real(8), intent(in)  :: temp\n")
f.write("real(8), intent(out) :: omegadot(ns) \n")
f.write('\n')
f.write('real(8) :: coi(ns), Tdiff \n')
f.write('real(8) :: M !< Third body\n')
f.write('integer :: is, T_i, Tint(2)\n')
f.write('real(8) :: prodf(1:'+str(nr)+'), prodb(1:'+str(nr)+')\n')
f.write('real(8) :: k(2) !< Falloff rate coefficients\n')
f.write('\n')
f.write('\n')

# Writing preliminary ops
f.write('do is = 1, ns \n')
f.write(' coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3\n')
f.write('enddo \n')
    
f.write('T_i = int(temp) \n')
f.write('Tdiff  = temp-T_i \n')
f.write('Tint(1) = T_i \n')
f.write('Tint(2) = T_i + 1 \n')

# Writing of prodf and prodb
ir=0; ir_arrhenius=0; ir_troe=0; ir_lindemann=0
for R in all_reactions:
  f.write('! reac n. '+str(ir+1)+': '+ str(R.equation)+'\n')
    
  # Standard and three-body Arrhenius
  if (R.reaction_type=='three-body-Arrhenius' or R.reaction_type=='Arrhenius'):
    ir_arrhenius += 1
    string_reac='prodf('+str(ir+1)+')=f_kf('+str(ir_arrhenius)+',Tint,Tdiff)'
    string_prod='prodb('+str(ir+1)+')=f_kb('+str(ir_arrhenius)+',Tint,Tdiff)'

    for react in R.reactants:
        ireact = find_indices(species, lambda e: e == react)[0]
        ni = gas.reactant_stoich_coeff(react, ir)

        if ni > 0:
            var = f'coi({ireact+1})'
            mult_expr = stoich_to_mult(var, ni)
            if mult_expr != '1.0':
                string_reac += f'*({mult_expr})'

    for prod in R.products:
        iprod = find_indices(species, lambda e: e == prod)[0]
        ni = gas.product_stoich_coeff(prod, ir)

        if ni > 0:
            var = f'coi({iprod+1})'
            mult_expr = stoich_to_mult(var, ni)
            if mult_expr != '1.0':
                string_prod += f'*({mult_expr})'
    
    if R.reaction_type == 'three-body-Arrhenius':
      string_reac=string_reac+'*M'
      string_prod=string_prod+'*M'
      # by default all species have efficiency 1.0
      efficiencies = R.third_body.efficiencies
      use_eff = len(efficiencies) > 0

      if use_eff:
          M_terms = []
          for isp, sp in enumerate(species):
            if sp in efficiencies:
              eff = efficiencies[sp]
            else:
              eff = R.third_body.efficiency(sp)
            coeff = format_eff(eff)
            if coeff != '*0': M_terms.append(f'coi({isp+1}){coeff}')
          Mstring = 'M=' + '+'.join(M_terms)
      else:
          Mstring = f'M=sum(coi(1:{len(species)}))'

      f.write(Mstring+'\n')

  elif (R.reaction_type=='falloff-Troe'):
    ir_troe += 1
    # by default all species have efficiency 1.0
    efficiencies = R.third_body.efficiencies
    use_eff = len(efficiencies) > 0

    if use_eff:
      M_terms = []
      for isp, sp in enumerate(species):
        eff = efficiencies.get(sp, 1.0)
        coeff = format_eff(eff)
        M_terms.append(f'coi({isp+1}){coeff}')
      Mstring = 'M=' + '+'.join(M_terms)
    else:
      Mstring = f'M=sum(coi(1:{len(species)}))'

    f.write(Mstring+'\n')
    f.write(f"k = f_k_troe({ir_troe},Tint,Tdiff,M)\n")
    string_reac='prodf('+str(ir+1)+')=k(1)'
    string_prod='prodb('+str(ir+1)+')=k(2)'

    for react in R.reactants:
        ireact = find_indices(species, lambda e: e == react)[0]
        ni = gas.reactant_stoich_coeff(react, ir)

        if ni > 0:
            var = f'coi({ireact+1})'
            mult_expr = stoich_to_mult(var, ni)
            if mult_expr != '1.0':
                string_reac += f'*({mult_expr})'

    for prod in R.products:
        iprod = find_indices(species, lambda e: e == prod)[0]
        ni = gas.product_stoich_coeff(prod, ir)

        if ni > 0:
            var = f'coi({iprod+1})'
            mult_expr = stoich_to_mult(var, ni)
            if mult_expr != '1.0':
                string_prod += f'*({mult_expr})'

  elif (R.reaction_type=='falloff-Lindemann'):
    ir_lindemann += 1
    # by default all species have efficiency 1.0
    efficiencies = R.third_body.efficiencies
    use_eff = len(efficiencies) > 0

    if use_eff:
      M_terms = []
      for isp, sp in enumerate(species):
        eff = efficiencies.get(sp, 1.0)
        coeff = format_eff(eff)
        M_terms.append(f'coi({isp+1}){coeff}')
      Mstring = 'M=' + '+'.join(M_terms)
    else:
      Mstring = f'M=sum(coi(1:{len(species)}))'

    f.write(Mstring+'\n')
    f.write(f"k = f_k_lindemann({ir_lindemann},Tint,Tdiff,M)\n")
    string_reac='prodf('+str(ir+1)+')=k(1)'
    string_prod='prodb('+str(ir+1)+')=k(2)'

    for react in R.reactants:
        ireact = find_indices(species, lambda e: e == react)[0]
        ni = gas.reactant_stoich_coeff(react, ir)

        if ni > 0:
            var = f'coi({ireact+1})'
            mult_expr = stoich_to_mult(var, ni)
            if mult_expr != '1.0':
                string_reac += f'*({mult_expr})'

    for prod in R.products:
        iprod = find_indices(species, lambda e: e == prod)[0]
        ni = gas.product_stoich_coeff(prod, ir)

        if ni > 0:
            var = f'coi({iprod+1})'
            mult_expr = stoich_to_mult(var, ni)
            if mult_expr != '1.0':
                string_prod += f'*({mult_expr})'

  f.write(string_reac+'\n')
  f.write(string_prod+'\n')
  
  ir=ir+1

# Writing of omegadot
f.write('! species source terms\n')
isp=1
for spec in species:
  # check if species is present in at least one reaction
  present=False 
  for R in all_reactions:
    if (spec in R.reactants) or (spec in R.products):
      present=True
       
  if (present):
    string='omegadot('+str(isp)+')=Wm_tab('+str(isp)+')*('
    ir=0
    for R in all_reactions:
      if (spec in R.reactants) or (spec in R.products):
        nir = gas.reactant_stoich_coeff(spec, ir)
        nip = gas.product_stoich_coeff(spec, ir)
        string=string+'+('+str(nip)+'-'+str(nir)+')*(prodf('+str(ir+1)+')-prodb('+str(ir+1)+'))'

      ir=ir+1
    string=string+(')')
  else:
    string='omegadot('+str(isp)+')=0.d0'
    
  f.write(string+'\n')
  
  isp=isp+1
  
f.write('end subroutine '+mechname+'\n')
f.write('end module '+mechname+'_mod\n')