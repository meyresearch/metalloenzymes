import multiprocessing.pool

def test_funct(x1):
    print(x1)



with multiprocessing.pool.Pool() as pool:
    pool.map(test_funct, ["lig1", "lig2"])
