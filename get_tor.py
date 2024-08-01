#!/usr/bin/env python3
#
# This scripts works only with the prmtop format introduced in Amber7

from sys import argv, exit

# DEFAULT VALUES
get_heavy=False  # get only those torsions involving heavy atoms
get_sidechain=False  # get only those torsions that are not backbone
get_backbone=False  # get only phi, psi, and omega angles
get_puck=False  # use puckering angles for 5-membered rings in PRO and HYP residues
noMet=False  # Exclude torsions for Methyl/Ammonium groups
select=False  # Select only the angles of a given set of residues

# READING OPTIONS
c = 0  # counter
if len(argv[0:]) == 1:
    print('Usage: get_tor.py [OPTIONS] amber_topology.top')
    print('for more info use the option -help')
    exit()
for i in range(1, len(argv[0:])):
  if argv[i] == '-help':
        print('\nUsage of get_tor.py:')
        print('\nSINOPSIS:')
        print('           get_tor.py [OPTIONS] amber_topology.top [> input_of_ptraj]')
        print('\nOPTIONS:')
        print('            -h      Select only torsion angles involving only heavy atoms')
        print('                    (i.e., not only Hydrogen)')
        print('\n          -b      Select only torsion angles defining the backbone')
        print('                    conformation of polypeptide/DNA molecules')
        print('                    (Note that -h -b is partially redundant)')
        print('\n          -s        Opposite to -b, select only the torsions that are not involved')
        print('                    in backbone conformation (i.e., side-chain torsions)')
        print('                    (Note that -b and -s are incompatible options)')
        print('\n          -sel     SELECTION')
        print('                    Select only torsions involving a given set of residues, for example:')
        print('                    -sel 5-20        selects torsions in residues from 5 to 20')
        print('                    -sel 1,5-20,33   as above, but including two more residues: 1 and 33')
        print('\n          -puck   The conformational state of 5-membered rings in PRO and HYP residues')
        print('                    is described by puckering angles (Pople-Cremer convention).')
        print('                    For DNA sugars, multipucker cpptraj method is selected.')
        print('\n          -noMet    Torsional angles involving Methyl/Ammonium groups are not considered')
        print('\n          -help     Print this quick help.')
        print('\nEXAMPLES:')
        print('>> get_tor.py peptide6.top > input_torsions.ptraj')
        print('>> get_tor.py -h -s -sel 3,4,7-10  peptide4.top > input.ptraj')
        exit()
  if argv[i]=='-sel':
     Selection_arg=argv[i+1]
     select=True
     c=c+2
  if argv[i]=='-b' :
     get_backbone=True
     c=c+1
  if argv[i]=='-s' :
     get_sidechain=True
     c=c+1
  if argv[i]=='-h' :
     get_heavy=True
     c=c+1
  if argv[i]=='-puck' :
     get_puck=True
     c=c+1
  if argv[i]=='-noMet' :
     noMet=True
     c=c+1
  if i==(len(argv[0:])-1):
     topology=open(argv[i],"r")
     c=c+1

# CHECKING OPTIONS
if (len(argv[0:])-1)!=c:
  print('\nERROR1: Check the options')
  exit()
if get_backbone==True and get_sidechain==True:
  print('\nERROR2: Incompatible options (-s) and (-b)')
  print('If you want to select both sidechain and backbone, do not specify anything')
  exit()

print ('trajin PUT-YOUR-TRAYECTORY-HERE')

#MAKING THE LIST OF SELECTED RESIDUES
temp_list1=[]
temp_list2=[]
SelectedResidues=[]
if (select):
    temp_list1=Selection_arg.split(',')
    for i in range(len(temp_list1)):
        temp_list2=temp_list1[i].split('-')
        if (len(temp_list2)>2):
            print ("Bad Argument",temp_list2)
            exit()
        elif (len(temp_list2)==2):
            for j in range( int(temp_list2[0]),(int(temp_list2[1])+1) ):
                SelectedResidues.append(j)
        else:
            SelectedResidues.append(int(temp_list1[i]))
        

# READING THE TOPOLOGY
begin_dih_noH=False
begin_dih_withH=False
begin_atom_name=False
begin_res_label=False
begin_residue_p=False
begin_bnd_H=False
begin_atnum=False
atom_name=[]
res_label=[]
temp_list=[]
dih_list=[]
bnd_H=[]
ibnd_H=[]
jbnd_H=[]
conn_H=[]
atnum=[]
res_pointer=[]
#Reading the atom names
for line in topology:
    temp_list=line.split()
    if line.find('ATOM_NAME') != -1:
        begin_atom_name=True
    if begin_atom_name==True and (line[0]!='%'):
        for i in range(len(temp_list)):
          while len(temp_list[i]) > 4:
            atom_name.append(temp_list[i][0:4])
            temp_list[i]=temp_list[i][4:]
          atom_name.append(temp_list[i])
    if line.find('FLAG CHARGE') !=-1:
       break

# Reading atomic numbers
for line in topology:
	temp_list=line.split()
	if line.find('ATOMIC_NUMBER') !=-1:
		begin_atnum=True
	if begin_atnum==True and (line[0]!='%'):
		for i in range(len(temp_list)):
			atnum.append(int(temp_list[i]))
	if line.find('MASS') !=-1:
		break
natom=len(atnum)
# print ('# Natom = ',natom)
#Reading the residue label 
for line in topology:
    temp_list=line.split()
    if line.find('FLAG RESIDUE_LABEL') != -1:
       begin_res_label=True
    if begin_res_label==True and (line[0]!='%'):
        for i in range(len(temp_list)):
          while len(temp_list[i]) > 4:
            res_label.append(temp_list[i][0:4])
            temp_list[i]=temp_list[i][4:]
          res_label.append(temp_list[i])
    if line.find('FLAG RESIDUE_POINTER') !=-1:
       topology.seek(0)
       break

#Reading the residue pointer
for line in topology:
    temp_list=line.split()
    if line.find('FLAG RESIDUE_POINTER') != -1:
        begin_residue_p=True
    if begin_residue_p and (line[0]!='%'):
        for i in range(len(temp_list)):
          res_pointer.append(int(temp_list[i]))
    if line.find('BOND_FORCE_CONSTANT') !=-1:
        res_pointer.append(len(atom_name)+1)
        break

#Reading bonds including hydrogen
for line in topology:
    temp_list=line.split()
    if line.find('BONDS_INC_HYDROGEN') !=-1:
         begin_bnd_H=True
    if begin_bnd_H==True and (line[0]!='%'):
        for i in range(len(temp_list)):
          bnd_H.append(int(temp_list[i]))
    if line.find('BONDS_WITHOUT_HYDROGEN') !=-1:
       break
nbnd_H= int (len(bnd_H)/3)
# print ('# Total number of Bonds including Hydrogen',nbnd_H)
j=0
for i in range(nbnd_H):
	ibnd_H.append(int( bnd_H[j]/3)+1  )
	jbnd_H.append( int( bnd_H[j+1]/3 ) +1)
	j=j+3

for i in range(natom):
	conn_H.append(0)

for i in range(nbnd_H):
	conn_H[ibnd_H[i]-1]=conn_H[ibnd_H[i]-1]+1
	conn_H[jbnd_H[i]-1]=conn_H[jbnd_H[i]-1]+1

#for i in range(natom):
#if ( atnum[i] == 6 ) and (conn_H[i] >= 3 ):
#	print( '# Methyl group at carbon ', i+1(

#Reading torsions without hydrogen
for line in topology:
    if line.find('DIHEDRALS_WITHOUT_HYDROGEN') != -1:
         begin_dih_noH=True
    if begin_dih_noH==True and (line[0]!='%'):
       if len(line.split())==10 and (int(line.split()[3])>=0):
           dih_list.append( line.split()[0:4] )
       if len(line.split())==10 and (int(line.split()[8])>=0):
           dih_list.append(line.split()[5:9])
       if len(line.split())==5 and (int(line.split()[3])>=0):
           dih_list.append(line.split()[0:4])
    if line.find('EXCLUDED_ATOMS_LIST') !=-1:
       break

#Reading torsions including hydrogens
topology.seek(0)
if (get_heavy==False):
 for line in topology:
    if line.find('DIHEDRALS_INC_HYDROGEN') != -1:
         begin_dih_withH=True
    if begin_dih_withH==True and (line[0]!='%'):
       if len(line.split())==10 and (int(line.split()[3])>=0):
           dih_list.append(line.split()[0:4])
       if len(line.split())==10 and (int(line.split()[8])>=0):
           dih_list.append(line.split()[5:9])
       if len(line.split())==5 and (int(line.split()[3])>=0):
           dih_list.append(line.split()[0:4])
    if line.find('DIHEDRALS_WITHOUT_HYDROGEN') !=-1:
       break

topology.close()

# Converting str into int
for i in range(len(dih_list)):
    for j in range(len(dih_list[i])):
        tmp_dih= int( dih_list[i][j] )
        dih_list[i][j]=tmp_dih


# GETTING RESNUMER[ATOMNUMBER]
ResNum=[]
for i in range(len(res_pointer)-1):
	for j in range(res_pointer[i],res_pointer[i+1]):
           ResNum.append(i+1)
SelectedAtoms=[]
for i in SelectedResidues :
    for j in range(len(atom_name)):
        if ResNum[j]==i:
            SelectedAtoms.append(j+1)  #Making the dih_list of selected atoms
# print('SelectedAtoms=',SelectedAtoms)

# MAKING THE FIRST LIST  OF TORSIONS
for i in range(len(dih_list)):           #Making the first dih_list of torsions
    dih_list[i][0]=abs(int(dih_list[i][0]/3))+1
    dih_list[i][1]=abs(int(dih_list[i][1]/3))+1
    dih_list[i][2]=abs(int(dih_list[i][2]/3))+1
    dih_list[i][3]=abs(int(dih_list[i][3]/3))+1

# OBTAINING A NEW LIST WITHOUT REDUNDANCY
new_dih_list=[]
k=0
for i in range(len(dih_list)):
    belong=False
#   Polypeptide
    belongs2back=False
    if     (atom_name[dih_list[i][1]-1]=='C' and atom_name[dih_list[i][2]-1]=='N') or \
           (atom_name[dih_list[i][2]-1]=='C' and atom_name[dih_list[i][1]-1]=='N') or \
           (atom_name[dih_list[i][1]-1]=='N' and atom_name[dih_list[i][2]-1]=='CA') or \
           (atom_name[dih_list[i][2]-1]=='N' and atom_name[dih_list[i][1]-1]=='CA') or \
           (atom_name[dih_list[i][1]-1]=='C' and atom_name[dih_list[i][2]-1]=='CA') or \
           (atom_name[dih_list[i][2]-1]=='C' and atom_name[dih_list[i][1]-1]=='CA'):
           belongs2back=True
#   DNA
    if     (atom_name[dih_list[i][1]-1]=='P' and atom_name[dih_list[i][2]-1]=='O5\'') or \
           (atom_name[dih_list[i][2]-1]=='O5\'' and atom_name[dih_list[i][1]-1]=='P') or \
           (atom_name[dih_list[i][1]-1]=='O5\'' and atom_name[dih_list[i][2]-1]=='C5\'') or \
           (atom_name[dih_list[i][2]-1]=='C5\'' and atom_name[dih_list[i][1]-1]=='O5\'') or \
           (atom_name[dih_list[i][1]-1]=='C5\'' and atom_name[dih_list[i][2]-1]=='C4\'') or \
           (atom_name[dih_list[i][2]-1]=='C4\'' and atom_name[dih_list[i][1]-1]=='C5\'') or \
           (atom_name[dih_list[i][1]-1]=='C4\'' and atom_name[dih_list[i][2]-1]=='C3\'') or \
           (atom_name[dih_list[i][2]-1]=='C3\'' and atom_name[dih_list[i][1]-1]=='C4\'') or \
           (atom_name[dih_list[i][1]-1]=='C3\'' and atom_name[dih_list[i][2]-1]=='O3\'') or \
           (atom_name[dih_list[i][2]-1]=='O3\'' and atom_name[dih_list[i][1]-1]=='C3\'') or \
           (atom_name[dih_list[i][1]-1]=='O3\'' and atom_name[dih_list[i][2]-1]=='P') or \
           (atom_name[dih_list[i][2]-1]=='P' and atom_name[dih_list[i][1]-1]=='O3\''):
           belongs2back=True

    if k>=0:
        for j in range(len(new_dih_list)): #Looking for redundancy in dih_list
            if (dih_list[i][1]==new_dih_list[j][1] and dih_list[i][2]==new_dih_list[j][2]) or \
            (dih_list[i][1]==new_dih_list[j][2] and dih_list[i][2]==new_dih_list[j][1]):
                 belong=True
        if (not belong):                    #Obtaining a new dih_list without redundancy
              if (get_backbone) and belongs2back :
                  new_dih_list.append(dih_list[i])
                  k+=1
              elif (get_sidechain) and (not belongs2back):
                 new_dih_list.append(dih_list[i])
                 k+=1
              elif (not get_backbone) and (not get_sidechain):
                 new_dih_list.append(dih_list[i])
                 k+=1

# FILTERING OUT TORSIONS INVOLVING METHYL GROUPS
if noMet:
    nxh3=0	
    new_dih_list2=[]
    for i in range(len(new_dih_list)):
        belongs1=True 
        belongs2=True
        iat=new_dih_list[i][1] 
        jat=new_dih_list[i][2] 
        if ( atnum[iat-1] == 6 ) and (conn_H[iat-1] >= 3 ):
            belongs1=False 
        if ( atnum[jat-1] == 6 ) and (conn_H[jat-1] >= 3 ):
            belongs2=False
        if ( atnum[iat-1] == 7 ) and (conn_H[iat-1] >= 3 ):
            belongs1=False 
        if ( atnum[jat-1] == 7 ) and (conn_H[jat-1] >= 3 ):
            belongs2=False
        if belongs1 and belongs2:
            new_dih_list2.append(new_dih_list[i])
        else:
            nxh3=nxh3+1

    new_dih_list=[]
    new_dih_list=new_dih_list2   
    print( '# Filtered ',nxh3,' torsions involving -CH3/-NH3+ groups')

# PICKING UP THE SELECTED ATOMS WITH THE OPTION "-sel"
if select:
    new_dih_list2=[]
    for i in range(len(new_dih_list)):
        belongs1=False
        belongs2=False
        for j in SelectedAtoms:
            if new_dih_list[i][1]==j:
                belongs1=True
            if new_dih_list[i][2]==j:
                belongs2=True
        if belongs1 and belongs2:
            new_dih_list2.append(new_dih_list[i])
    new_dih_list=[]   
    new_dih_list=new_dih_list2   
#   print('new_dih_list=',new_dih_list2)

#FILTERING OUT CHI DIHEDRAL ANGLES IN HYP/PRO RINGS
if (get_puck) :
   new_dih_list2=[]
   for i in range(len(new_dih_list)):
         belongs1=True
         atom_1=atom_name[new_dih_list[i][1]-1] 
         atom_2=atom_name[new_dih_list[i][2]-1] 
         res_1=res_label[ResNum[new_dih_list[i][1]-1]-1]
         res_2=res_label[ResNum[new_dih_list[i][2]-1]-1]
         if atom_1 == 'CB' and atom_2 == 'CA'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if atom_1 == 'CG' and atom_2 == 'CB'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if atom_1 == 'CD' and atom_2 == 'CG'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if atom_1 == 'N' and atom_2 == 'CD'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if (belongs1) : new_dih_list2.append(new_dih_list[i])
   new_dih_list=[]   
   new_dih_list=new_dih_list2   

if (get_puck) :
   new_dih_list2=[]
   for i in range(len(new_dih_list)):
         belongs1=True
         atom_1=atom_name[new_dih_list[i][1]-1] 
         atom_2=atom_name[new_dih_list[i][2]-1] 
         res_1=res_label[ResNum[new_dih_list[i][1]-1]-1]
         res_2=res_label[ResNum[new_dih_list[i][2]-1]-1]
         if atom_1 == 'CA' and atom_2 == 'CB'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if atom_1 == 'CB' and atom_2 == 'CG'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if atom_1 == 'CG' and atom_2 == 'CD'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if atom_1 == 'CD' and atom_2 == 'N'  and \
            (res_1 == 'PRO' or res_1 == 'HYP') and \
            (res_2 == 'PRO' or res_2 == 'HYP') : belongs1=False
         if (belongs1) : new_dih_list2.append(new_dih_list[i])
   new_dih_list=[]   
   new_dih_list=new_dih_list2   

#FILTERING OUT NU DIHEDRAL ANGLES IN NUCLEIC ACID SUGARS
if (get_puck) :
   new_dih_list2=[]
   for i in range(len(new_dih_list)):
         belongs1=True
         atom_1=atom_name[new_dih_list[i][1]-1] 
         atom_2=atom_name[new_dih_list[i][2]-1] 
         if ( atom_1 == 'O4\'' and atom_2 == 'C1\''  ) : belongs1=False
         if ( atom_1 == 'C1\'' and atom_2 == 'O4\''  ) : belongs1=False
         if ( atom_1 == 'C1\'' and atom_2 == 'C2\''  ) : belongs1=False
         if ( atom_1 == 'C2\'' and atom_2 == 'C1\''  ) : belongs1=False
         if ( atom_1 == 'C2\'' and atom_2 == 'C3\''  ) : belongs1=False
         if ( atom_1 == 'C3\'' and atom_2 == 'C2\''  ) : belongs1=False
         if ( atom_1 == 'C3\'' and atom_2 == 'C4\''  ) : belongs1=False
         if ( atom_1 == 'C4\'' and atom_2 == 'C3\''  ) : belongs1=False
         if ( atom_1 == 'C4\'' and atom_2 == 'O4\''  ) : belongs1=False
         if ( atom_1 == 'O4\'' and atom_2 == 'C4\''  ) : belongs1=False
         if (belongs1) : new_dih_list2.append(new_dih_list[i])
   new_dih_list=[]   
   new_dih_list=new_dih_list2   
	        

# PRINTING OUT TORSION DATA
ntor=len(new_dih_list)
print ('# Total number of torsion angles ',ntor)
for i in range(len(new_dih_list)):
    if (i+1)<10:                    name='d'+'000'+str(i+1)+'.dat'
    if (i+1)>=10 and (i+1)<100:     name='d'+'00'+str(i+1)+'.dat'
    if (i+1)>=100 and (i+1)<1000:   name='d'+'0'+str(i+1)+'.dat'
    if (i+1)>=1000:                 name='d'+str(i+1)+'.dat'
    alias1=''
    alias2=''
    for k in [0,1,2,3] : 
           alias1=alias1+res_label[ResNum[new_dih_list[i][k]-1]-1]+'-'+str(ResNum[new_dih_list[i][k]-1])+'@'+atom_name[new_dih_list[i][k]-1]
           alias2=alias2+res_label[ResNum[new_dih_list[i][k]-1]-1]+'_'+str(ResNum[new_dih_list[i][k]-1])+'_'+atom_name[new_dih_list[i][k]-1]
           if (k < 3 ) : alias2=alias2+'-'
    print( '# Atoms involved in torsion ',name[0:5],' ',alias1)
    print ('dihedral ',alias2,'@'+str(new_dih_list[i][0]),'@'+str(new_dih_list[i][1]),'@'+str(new_dih_list[i][2]),'@'+str(new_dih_list[i][3]), ' out ',name)

# PRINTING OUT PUCKERING DATA HYP/PRO
npuck=0
if (get_puck) :
   for i in range(len(res_label)):
         if (res_label[i] == 'PRO') or (res_label[i] == 'HYP' ) : npuck=npuck+1
   print('# Total number of PRO/HYP puckering angles ',npuck)
   ipuck=ntor
   for i in range(len(res_label)):
         if (res_label[i] == 'PRO') or (res_label[i] == 'HYP' )  : 
              ipuck=ipuck+1
              if (ipuck)<10:                       name='d'+'000'+str(ipuck)+'.dat'
              if (ipuck) >=10 and (ipuck)<100:     name='d'+'00'+str(ipuck)+'.dat'
              if (ipuck)>=100 and (ipuck)<1000:    name='d'+'0'+str(ipuck)+'.dat'
              if (ipuck)>=1000:                    name='d'+str(ipuck)+'.dat'
              print( '# Atoms involved in puckering ',name[0:5],' ',res_label[i]+'-'+str(i+1)+'@N,CA,CB,CG,CD')
              alias=res_label[i]+'_'+str(i+1)+'_ring_pucker_angle'
              print( 'pucker ',alias,' ',':'+str(i+1)+'@N',':'+str(i+1)+'@CA',':'+str(i+1)+'@CB',':'+str(i+1)+'@CG',':'+str(i+1)+'@CD', ' out',name,' cremer')

# PRINTING OUT PUCKERING DATA NUCLEIC ACIDS 
npuck_nuc=0
if (get_puck) :
   for i in range(len(res_label)):
         if ( 'DT' in res_label[i] ) or ( 'DA' in res_label[i] ) \
         or ( 'DG' in res_label[i] ) or ( 'DC' in res_label[i] ) : npuck_nuc=npuck_nuc+1
   print('# Total number of nucleic puckering angles ',npuck_nuc)
   ipuck=ntor+npuck
   for i in range(len(res_label)):
         if ( 'DT' in res_label[i] ) or ( 'DA' in res_label[i] ) \
         or ( 'DG' in res_label[i] ) or ( 'DC' in res_label[i] ) : 
              ipuck=ipuck+1
              if (ipuck)<10:                       name='d'+'000'+str(ipuck)+'.dat'
              if (ipuck) >=10 and (ipuck)<100:     name='d'+'00'+str(ipuck)+'.dat'
              if (ipuck)>=100 and (ipuck)<1000:    name='d'+'0'+str(ipuck)+'.dat'
              if (ipuck)>=1000:                    name='d'+str(ipuck)+'.dat'
              print( '# Atoms involved in puckering ',name[0:5],' ',res_label[i]+'-'+str(i+1)+'@O4\',C1\',C2\',C3\',C4\'')
              alias=res_label[i]+'_'+str(i+1)+'_ring_pucker_angle'
              print( 'multipucker ',alias,' nucleic resrange ',str(i+1)+'-'+str(i+1),' out ',name)

# ARRANGING THE UNIQUE LIST OF CENTRAL ATOMS
central_dih_list=[] 
for i in range(len(new_dih_list)):
	central_dih_list.append(new_dih_list[i][1])

# for puckering angles we choose CB/C1' and CD/C4' as central atoms 
if ( get_puck) :
    for i in range(len(res_label)):
        if (res_label[i] == 'PRO') or (res_label[i] == 'HYP' ) :
            for j in range(res_pointer[i],res_pointer[i+1]):
                if atom_name[j-1] == 'CB' :  central_dih_list.append(j)

    for i in range(len(res_label)):
        if ( 'DT' in res_label[i] ) or ( 'DA' in res_label[i] ) \
        or ( 'DG' in res_label[i] ) or ( 'DC' in res_label[i] ) : 
            for j in range(res_pointer[i],res_pointer[i+1]):
                if atom_name[j-1] == 'C1\'' :  central_dih_list.append(j)

for i in range(len(new_dih_list)):
	central_dih_list.append(new_dih_list[i][2])

# for puckering angles we choose CB/C1' and CD/C4' as central atoms 
if ( get_puck) :
    for i in range(len(res_label)):
        if (res_label[i] == 'PRO') or (res_label[i] == 'HYP' ) :
            for j in range(res_pointer[i],res_pointer[i+1]):
                if atom_name[j-1] == 'CD' :  central_dih_list.append(j)

    for i in range(len(res_label)):
        if ( 'DT' in res_label[i] ) or ( 'DA' in res_label[i] ) \
        or ( 'DG' in res_label[i] ) or ( 'DC' in res_label[i] ) : 
            for j in range(res_pointer[i],res_pointer[i+1]):
                if atom_name[j-1] == 'C4\'' :  central_dih_list.append(j)

uniq_dih_list=[] 
uniq_dih_list.append(central_dih_list[0])
for i in range(len(central_dih_list))[1:len(central_dih_list)]:
    check=1
    for j in range(len(uniq_dih_list)):
    	if central_dih_list[i] == uniq_dih_list[j]:
            check=0
    if (check): 
        uniq_dih_list.append(central_dih_list[i])

uniq_dih_list.sort()

# BUILDING THE REDUCED LIST OF CENTRAL ATOMS 
id_central_dih_list=[]
for i in range(len(central_dih_list)):
	k=0
	for j in range(len(uniq_dih_list)):
		if central_dih_list[i] == uniq_dih_list[j]: k=j
	id_central_dih_list.append(k)

# WRITING THE CENTRAL ATOMS OF TORSION/PUCKERING ANGLES USING THE REDUCED LIST 
info=open("atoms_in_tor.info","w")
for i in range(ntor+npuck+npuck_nuc):
    info.write(str(id_central_dih_list[i]+1)+'  '+str(id_central_dih_list[i+ntor+npuck+npuck_nuc]+1)+'\n')

print ('matrix dist \\')
for i in range(len(uniq_dih_list)):
	if (i==0) :
		print  ('@'+str(uniq_dih_list[i]+1)+',\\')
	elif (i<len(uniq_dih_list)-1):
		print ( str(uniq_dih_list[i]+1)+',\\')
	else:
		print ( str(uniq_dih_list[i]+1)+' \\')
for i in range(len(uniq_dih_list)):
	if (i==0) :
		print ( '@'+str(uniq_dih_list[i]+1)+',\\')
	elif (i<len(uniq_dih_list)-1):
		print ( str(uniq_dih_list[i]+1)+',\\')
	else:
		print ( str(uniq_dih_list[i]+1)+' \\')
print (' out distance_matrix.dat' )

