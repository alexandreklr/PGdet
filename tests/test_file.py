from PGdet import atom_mapping, count_atoms, calculate_angle, pg_check_smiles
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
  smiles = ''
  for i in range(27):
    assert count_atoms(smiles) == i, 'Test failed : The wrong value was returned.'
    smiles += 'C'

def test_calculate_angle():
    assert calculate_angle([1,2,3],[4,5,6],[7,8,9]) == 180.0, 'Test failed : Wrong value was returned.'
    assert calculate_angle([3,-1,4],[2,2,2],[0,1,-3]) == 116.0499813323194, 'Test failed : Wrong value was returned.'
    assert calculate_angle([2,3,4],[5,6,7],[8,9,11]) == 171.9505330244716, 'Test failed : Wrong value was returned.'
    assert calculate_angle([0,0,0],[1,2,2],[2,1,3]) == 101.09580328313639, 'Test failed : Wrong value was returned.'
    assert calculate_angle([7,3,2],[6,4,1],[5,5,1]) == 144.73561031724535, 'Test failed : Wrong value was returned.'
    assert calculate_angle([2,3,4],[3,3,3],[2,2,2]) == 90.0, 'Test failed : Wrong value was returned.'
    assert calculate_angle([8,7,6],[5,5,5],[2,4,4]) == 165.23686737636845, 'Test failed : Wrong value was returned.'
  
def test_pg_check_smiles():
  try:
    test_smiles = ['CC', 'CCO', 'NN', 'c1ccccc1', 'CCCCCCCCCCCCCCCCCCCCCCCCC']
    for smile in test_smiles:
      pg_check_smiles(smile)
  except TypeError:
    assert False, 'Test failed : Error raised when not expected.'

  test_smiles = ['not_a_smiles', 'FS(F)(F)(F)(F)F', 'CCCCCCCCCCCCCCCCCCCCCCCCCC', '', 3, ('invalid' , 'set'), {'invalid' : 'dictionary'}, True, '[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Fe+2]']
  for smile in test_smiles:
    try :
      pg_check_smiles(smile)
      assert False, 'Test failed : Expected error not raised'
    except TypeError :
        continue
    
