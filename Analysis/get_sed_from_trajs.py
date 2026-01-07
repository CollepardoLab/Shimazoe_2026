import sys
import numpy as np
from calc_sed import calc_sed
# load in full traj files and compute the sed coeff of each frame


def get_frames(fname):
    frames=list()
    traj_file=open(fname)
    I=1
    frame=list()

    for line in traj_file:
        if line == "ITEM: TIMESTEP\n":
            print(I)
            if I>1:
                frame=np.array(frame)
                frames.append(frame)
                frame=list()
            I=I+1

        else: # try and read in line
            try:
                array_line = np.array(line.split(), dtype=float)
                if len(array_line) == 14:
                    frame.append(array_line)
            except:
                pass

    traj_file.close()
    # last frame
    frame=np.array(frame)
    frames.append(frame)
    return frames


for i in range(0,1):
    print(i)
    fname = '../equil_T0.dump'

    print("#"+fname)
    frames=get_frames(fname)
    ss=list()
    for frame in frames:
        s=calc_sed(frame[:,:])
        print(s)
        ss.append(s)

    with open('Sedimentation.txt', 'w') as f:
        f.write('{} {}\n'.format(np.mean(ss[-100:]), np.std(ss[-100:])))
    f.close()

    with open ('195_MobileH1.txt', 'w') as ff:
        for ele in ss[-100:]:
            ff.write('{}\n'.format(ele))
    ff.close()

     # Convert ss to numpy array
    ss = np.array(ss)
    # Calculate cumulative average
    cumulative_avg = np.cumsum(ss) / np.arange(1, len(ss) + 1)

    # Plot the data points and cumulative average
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(ss, label='Sedimentation Coefficient', marker='o')
    plt.plot(cumulative_avg, label='Cumulative Average', linestyle='--')
    plt.xlabel('Frame Index')
    plt.ylabel('Sedimentation Coefficient')
    plt.legend()
    plt.title('Sedimentation Coefficient vs Time')
    plt.grid()
    plt.show()
