with open('ligands_xyz.xyz', 'r') as fh:
    xyz_lines = fh.read().split('\n')

for i in range(len(xyz_lines)):

    line_split = xyz_lines[i].split(' ')
    if len(line_split) != 4:
        continue

    for j in range(1, len(line_split), 1):

        if 'e' in line_split[j]:
            line_split[j] = '{:6f}'.format(float(line_split[j]))

    xyz_lines[i] = ' '.join(line_split)

with open('ligands_xyz.xyz', 'w') as fh:
   fh.write('\n'.join(xyz_lines))
