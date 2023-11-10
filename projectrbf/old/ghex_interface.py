import numpy as np

from fixtures.context import *

import ghex.unstructured as unstructured

LEVELS=2 # number of levels vertically / need clarity and discussion

#all, outer, outer lids to be obtained from grid manager
#domains data structure to be replaced with the domain data
#subdomains = ranks = global owners of processes from grid manager

def mydomain_descriptor(mpi_cart_comm):
    comm = ghex.mpi_comm(mpi_cart_comm)
    ctx = ghex.context(comm, True)
    #assert ctx.size() == 4

    domain_desc = unstructured.domain_descriptor(
            ctx.rank(),
            domains[ctx.rank()]["all"],
            domains[ctx.rank()]["outer_lids"])

    # assert domain_desc.domain_id() == ctx.rank()
    # assert domain_desc.size() == len(domains[ctx.rank()]["all"])
    # assert domain_desc.inner_size() == len(domains[ctx.rank()]["inner"])

    #halo_gen = unstructured.halo_generator()
    halo_gen = unstructured.halo_generator_with_gids(domains[ctx.rank()]["outer"])

    pattern = unstructured.make_pattern(ctx, halo_gen, [domain_desc])

    co = unstructured.make_co(ctx) #communication object

    # make_feield function is creating the f vector in the Ac = f system
    def make_field():
        data = np.zeros([len(domains[ctx.rank()]["all"]), LEVELS], dtype=np.float64) #what is the levels dimension ?
        inner_set = set(domains[ctx.rank()]["inner"])
        all_list = domains[ctx.rank()]["all"]
        for x in range(len(all_list)):
            gid = all_list[x]
            for l in range(LEVELS):
                if gid in inner_set:
                    data[x, l] = ctx.rank()*1000 + 10*gid + l
                    #data[x,l] = sample data field
                else:
                    data[x, l] = -1

        field = unstructured.field_descriptor(domain_desc, data)
        return data, field


    d1, f1 = make_field()
    #d2, f2 = make_field()

    res = co.exchange([pattern(f1)])
    res.wait()

    return d1,f1



    # def check_field(data):
    #     inner_set = set(domains[ctx.rank()]["inner"])
    #     all_list = domains[ctx.rank()]["all"]
    #     for x in range(len(all_list)):
    #         gid = all_list[x]
    #         for l in range(LEVELS):
    #             if gid in inner_set:
    #                 assert data[x, l] == ctx.rank()*1000 + 10*gid + l
    #             else:
    #                 assert (data[x, l] - 1000*int((data[x, l])/1000)) == 10*gid + l
    #
    #     field = unstructured.field_descriptor(domain_desc, data)
    #     return data, field


    # check_field(d1)
    # check_field(d2)
