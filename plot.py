import matplotlib.pyplot as plt
x1 = []
y1 = []
y2 = []
with open('test_p_1', 'r') as f:
    line = f.readline()
    while line:
        s = line.split(' ')
        a = int(s[0])
        b = float(s[1])
        c = float(s[2])
        x1.append(s)
        y1.append(b)
        y2.append(c)
        line = f.readline()
f.closed
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x1, y1, 'b-')
ax1.set_xlabel('atom')
ax1.set_ylabel('ccij', color='b')
for t1 in ax1.get_yticklabels():
    t1.set_color('b')

ax2 = ax1.twinx()
ax2.plot(x1, y2, 'r-')
ax2.set_ylabel('k-l', color='r')
for t1 in ax2.get_yticklabels():
    t1.set_color('r')

plt.show()
