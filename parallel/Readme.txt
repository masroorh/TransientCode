***
*** using million elements mesh ***
***
change chunksize from 64 to 4096 in 3d-darcy.f
change isize=1200000

***
*** running transient algorithm
***
change trans from 0 to 1 in 3d-darcy.f

***
*** running metis programs
***
        gfortran 3d-darcy-metis.f 3d-darcy-elmt.f -o 3ddarcy-metis -w -fopenmp

        /* mesh to graph conversion */
        m2gmetis -gtype=nodal geom_Darcy_xxxxx_elem.dat geom_Darcy_xxxxx_graph.dat

        /* multilevel mesh reordering */
        ndmetis geom_Darcy_xxxxx_graph.dat

        cp geom_Darcy_xxxxx_graph.dat.iperm metis_geom_Darcy.dat

        cp geom_Darcy_xxxxx.dat geom_Darcy.dat

        ./3ddarcy-metis

        /* multhreaded multilevel mesh reordering */
        time mtmetis -pnd geom_Darcy_20536.graph geom_Darcy_20536.perm

