import multiprocessing
import threading
from sys import stdout
from collections import Counter

#http://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class/21345308
#https://docs.python.org/3/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes

def listener_builder(sequence_dict, outdir):

    def listener_process(q_out):
        data = Counter()
        while True:
            result = q_out.get()
            if result is None:
                break
            else:
                data[result]+=1

        for offset in ['P0', 'P1', 'P2', 'AP0', 'AP1', 'AP2']:
            orient = offset[:-1]
            frame = int(offset[-1])
            for gene, sequence in sequence_dict.items():
                keys = range(1, int((len(sequence)/3)))
                with open('%s/%s_%s_%d.count'%(outdir, gene, orient, frame), 'w') as f:
                    for k in keys:
                        count = data.get((gene, orient+'1', float(k), frame),0) + data.get((gene, orient+'2', float(k), frame), 0)
                        f.write("%d %d\n"%(k, count))
        return
    return listener_process

def fun(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put(f(x))

def parmap(f, X, listener_process, nprocs=multiprocessing.cpu_count()):
    q_in=multiprocessing.Queue(1000)
    q_out=multiprocessing.Queue(1000)

    proc=[multiprocessing.Process(target=fun, args=(f, q_in, q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon=True
        p.start()

    listener = threading.Thread(target=listener_process, args=(q_out,))
    listener.daemon=True
    listener.start()

    [q_in.put((i,x), timeout=3600) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    [p.join() for p in proc]
    q_out.put(None)

    listener.join()