import numpy as np
import math
import copy

def DoAutoGainControl(self, float timewin):
    """ Cython version of AGC algorithm.

        timewin: time half-window in seconds
        dt: sample interval in seconds
    """
    cdef float dt
    if self.rate is not None:
        dt = self.rate
    else:
        dt = 1e-8

    cdef list data = self.data.tolist()
    cdef int iwagc = int(round(timewin/dt))     # Half-window in samples
    cdef int nwin = iwagc*2 + 1
    cdef int i, j
    cdef float x, y, xk
    cdef float blocksum, scl_const
    cdef list X, Y
    cdef list trace = [0. for i in range(len(data))]
    cdef list agctrace = copy.copy(trace)

    for i in range(len(data[0])):

        if sum([abs(X[i]) for X in data]) != 0:
            nt = len(data)
            trace = map(lambda X:X[i], data)
            agctrace = [0. for x in trace]
            blocksum = sum([x**2 for x in trace[:iwagc]])

            # Compute the first half-window
            for j in range(iwagc):
                blocksum = blocksum + trace[j+iwagc]**2
                agctrace[j] = trace[j] / math.sqrt(blocksum / (j+iwagc))

            # Compute the middle, full-rms window
            for j in range(iwagc, nt-iwagc):
                blocksum = blocksum - trace[j-iwagc]**2
                blocksum = blocksum + trace[j+iwagc]**2
                try:
                    agctrace[j] = trace[j] / math.sqrt(blocksum / (2*iwagc))
                except ValueError:
                    # Sometimes there's a domain error here - should look into it
                    pass

            # Compute the last half-window
            for j in range(nt-iwagc, nt):
                blocksum = blocksum - trace[j-iwagc]**2
                try:
                    agctrace[j] = trace[j] / math.sqrt(blocksum / (nt-j+iwagc))
                except ValueError:
                    # Sometimes there's a domain error here - should look into it
                    pass

            # Scale to preserve average amplitude
            #scl_const = np.sqrt(np.mean(trace**2)) / np.sqrt(np.mean(agctrace**2))
            x = sum([xk**2 for xk in trace]) / float(len(trace))
            y = sum([xk**2 for xk in agctrace]) / float(len(agctrace))
            scl_const = math.sqrt(x) / math.sqrt(y)
            for j in range(len(agctrace)):
                agctrace[j] = agctrace[j] * scl_const

            # Copy array into new line
            self.data[:,i] = np.array(agctrace)

    return
