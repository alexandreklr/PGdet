from PGdet import atom_mapping, count_atoms, calculate_angle, pg_check_smiles, config_mol, molecule_angles, molecule_bonds

def test_atom_mapping():
  assert atom_mapping('CCO') == 'Cs', 'Test failed : The wrong point group was returned'
  assert atom_mapping('ICl') == 'Cinfv', 'Test failed : The wrong point group was returned'
  assert atom_mapping(['CC', 'C', 'c1ccccc1']) == ['D3d', 'Td', 'D6h'], 'Test failed : The wrong point group was returned'

  # Wrong input testing
  try :
    atom_mapping('wrong smiles')
    assert False, 'Test failed : Expected error not raised.'
  except TypeError:
    pass
  try :
    atom_mapping('FS(F)(F)(F)(F)F')
    assert False, 'Test failed : Expected error not raised.'
  except TypeError:
    pass

def test_count_atoms():
  pass

def test_calculate_angle():
  pass

def test_pg_check_smiles():
  pass

def test_config_mol():
  pass

def test_molecule_angles():
  pass

def test_molecules_bonds():
  
