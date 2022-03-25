
'''
# Find nearest crosslinking group (10, 11) to the broken bond
crosslinking = ['10','11']
checked_atoms = []
checked_bonds = []
broken = True
while broken:

    # check all atoms bonded to the first atom in uncrosslinked bond
    first_atom = bonds[uncrosslink]['atoms'][0]
    bonded_atoms = atoms[first_atom]['bonded']
    checked_atoms.append(first_atom)

    new_bond = False
    
    for a in bonded_atoms:
        if a not in checked_atoms:
            next_atom = atoms[a]

            if next_atom['type'] in crosslinking:
                print('Found a crosslinking group nearby!')
                broken = False
                break
            else:
                print('Checking next bonded atom...')
                checked_atoms.append(next_atom)

    if not broken:
        break

    # check all atoms bonded to the second atom in uncrosslinked bond
    second_atom = bonds[uncrosslink]['atoms'][1]
    bonded_atoms = atoms[second_atom]['bonded']
    checked_atoms.append(second_atom)
    
    for a in bonded_atoms:
        if a not in checked_atoms:
            next_atom = atoms[a]

            if next_atom['type'] in crosslinking:
                print('Found a crosslinking group nearby!')
                broken = False
                break
            else:
                print('Checking next bonded atom...')
                checked_atoms.append(next_atom)

    if not broken:
        break

    checked_bonds.append(uncrosslink)
    # Find a bond that has the first atom and has not been checked
    for bond in bonds:
        if first_atom in bonds[bond]['atoms'] and bond not in checked_bonds:
            uncrosslinked = bond
            new_bond = True
            break

    if not new_bond: # if no bond that has first atom and has not been checked 
        
        for bond in bonds: # look for a bond with second atom and has not been checked
            if second_atom in bonds[bond]['atoms'] and bond not in checked_bonds:
                uncrosslinked = bond
                new_bond = True
                break
'''



'''
for b in bonds:

    # Find nearest crosslinking group (10, 11) to the broken bond
    crosslinking = ['10','11']
    checked_atoms = []
    checked_bonds = []
    broken = bonds[b]['broken']
    uncrosslink = b
    while broken and not bonds[b]['delete']:

        # check all atoms bonded to the first atom in uncrosslinked bond
        print('Checking bond %s:' %(uncrosslink))
        print(bonds[uncrosslink])
        print()

        if bonds[uncrosslink]['atoms'][0] not in checked_atoms: # check to see if first atom has already been checked
            first_atom = bonds[uncrosslink]['atoms'][0]
            bonded_atoms = atoms[first_atom]['bonded']
            checked_atoms.append(first_atom)

            new_bond = False
            
            for a in bonded_atoms:
                if a not in checked_atoms:
                    next_atom = atoms[a]

                    if next_atom['type'] in crosslinking:
                        print('Found a crosslinking group nearby!')
                        print(first_atom,bonded_atoms)
                        print(a,atoms[a])
                        for bid in atoms[a]['bonds']: # find the crosslinking bond
                            if bonds[bid]['type'] == '5': 
                                bonds[bid]['broken'] = True
                                bonds[bid]['delete'] = True
                                break
                        print(bid,bonds[bid])
                        print()
                        broken = False
                        break
                    else:
                        checked_atoms.append(a)

            if not broken:
                break

            # check all atoms bonded to the second atom in uncrosslinked bond
            second_atom = bonds[uncrosslink]['atoms'][1]
            bonded_atoms = atoms[second_atom]['bonded']
            checked_atoms.append(second_atom)
            
            for a in bonded_atoms:
                if a not in checked_atoms:
                    next_atom = atoms[a]

                    if next_atom['type'] in crosslinking: # if the next atom is a crosslinking atom
                        print('Found a crosslinking group nearby!')
                        print(second_atom,bonded_atoms)
                        print(a,atoms[a])
                        for bid in atoms[a]['bonds']: # find the crosslinking bond
                            if bonds[bid]['type'] == '5': 
                                bonds[bid]['broken'] = True
                                bonds[bid]['delete'] = True
                                break
                        print(bid,bonds[bid])
                        print()
                        broken = False
                        break
                    else:
                        checked_atoms.append(a)

            if not broken:
                break

            checked_bonds.append(uncrosslink)
            # Find a bond that has the first atom and has not been checked
            for bond in bonds:
                if first_atom in bonds[bond]['atoms'] and bond not in checked_bonds:
                    uncrosslink = bond
                    new_bond = True
                    break

            if not new_bond: # if no bond that has first atom and has not been checked 
                
                for bond in bonds: # look for a bond with second atom and has not been checked
                    if second_atom in bonds[bond]['atoms'] and bond not in checked_bonds:
                        uncrosslink = bond
                        new_bond = True
                        break

        else: # if first atom has already been checked
            second_atom = bonds[uncrosslink]['atoms'][1]
            bonded_atoms = atoms[second_atom]['bonded']
            checked_atoms.append(second_atom)
            
            for a in bonded_atoms:
                if a not in checked_atoms:
                    next_atom = atoms[a]

                    if next_atom['type'] in crosslinking: # if the next atom is a crosslinking atom
                        print('Found a crosslinking group nearby!')
                        print(second_atom,bonded_atoms)
                        print(a,atoms[a])
                        for bid in atoms[a]['bonds']: # find the crosslinking bond
                            if bonds[bid]['type'] == '5': 
                                bonds[bid]['broken'] = True
                                bonds[bid]['delete'] = True
                                break
                        print(bid,bonds[bid])
                        print()
                        broken = False
                        break
                    else:
                        checked_atoms.append(a)

            if not broken:
                break

            checked_bonds.append(uncrosslink)

            for bond in bonds: # look for a bond with second atom and has not been checked
                if second_atom in bonds[bond]['atoms'] and bond not in checked_bonds:
                    uncrosslink = bond
                    new_bond = True
                    break
'''

'''
for b in bonds:

    # Find nearest crosslinking group (10, 11) to the broken bond
    crosslinking = ['10','11']
    checked_atoms = []
    checked_bonds = []
    broken = bonds[b]['broken']
    uncrosslink = b
    while broken and not bonds[b]['delete']:

        if uncrosslink not in checked_bonds: # if we have not checked this bond before

            # check all atoms bonded to the first atom in uncrosslinked bond
            print('Checking bond %s:' %(uncrosslink))
            print(bonds[uncrosslink])
            print()
            first_atom = bonds[uncrosslink]['atoms'][0]
            bonded_atoms = atoms[first_atom]['bonded']
            checked_atoms.append(first_atom)

            new_bond = False
            
            for a in bonded_atoms:
                if a not in checked_atoms:
                    next_atom = atoms[a]

                    if next_atom['type'] in crosslinking:
                        print('Found a crosslinking group nearby!')
                        print(first_atom,bonded_atoms)
                        print(a,atoms[a])
                        for bid in atoms[a]['bonds']: # find the crosslinking bond
                            if bonds[bid]['type'] == '5': 
                                bonds[bid]['broken'] = True
                                bonds[bid]['delete'] = True
                                break
                        print(bid,bonds[bid])
                        print()
                        broken = False
                        break
                    else:
                        checked_atoms.append(a)

            if not broken:
                break

            # check all atoms bonded to the second atom in uncrosslinked bond
            second_atom = bonds[uncrosslink]['atoms'][1]
            bonded_atoms = atoms[second_atom]['bonded']
            checked_atoms.append(second_atom)
            
            for a in bonded_atoms:
                if a not in checked_atoms:
                    next_atom = atoms[a]

                    if next_atom['type'] in crosslinking: # if the next atom is a crosslinking atom
                        print('Found a crosslinking group nearby!')
                        print(second_atom,bonded_atoms)
                        print(a,atoms[a])
                        for bid in atoms[a]['bonds']: # find the crosslinking bond
                            if bonds[bid]['type'] == '5': 
                                bonds[bid]['broken'] = True
                                bonds[bid]['delete'] = True
                                break
                        print(bid,bonds[bid])
                        print()
                        broken = False
                        break
                    else:
                        checked_atoms.append(a)

            if not broken:
                break

            checked_bonds.append(uncrosslink)
            # Find a bond that has the first atom and has not been checked
            for bond in bonds:
                if first_atom in bonds[bond]['atoms'] and bond not in checked_bonds:
                    previous_bond = uncrosslink
                    uncrosslink = bond
                    new_bond = True
                    break

            if not new_bond: # if no bond that has first atom and has not been checked 
                
                for bond in bonds: # look for a bond with second atom and has not been checked
                    if second_atom in bonds[bond]['atoms'] and bond not in checked_bonds:
                        previous_bond = uncrosslink
                        uncrosslink = bond
                        new_bond = True
                        break

        else: # if we have checked this bond before, need to pick a new one
            print('Already checked bond %s... Trying previous bond again...' %(uncrosslink))
            second_atom = bonds[previous_bond]['atoms'][1]
            for bid in atoms[second_atom]['bonds']:
                if bid not in checked_bonds:
                    uncrosslink = bid
                    checked_bonds = []
                    checked_atoms = []
                    break
'''


