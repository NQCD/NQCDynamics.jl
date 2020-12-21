
function H2_Cu_PES_out(input_path::String, coordinates_type::Int64=0, force_incl::Int64=1, constraint_num::Int64=0, generate_output::Bool=false)
    # EANN PES: if the code was crashed during last run, use the line below (otherwise Julia will most likely crush)
    #ccall((:deallocate_all_, "./nn.dylib"),Cvoid ,())
    # EANN PES: call init_pes function

    #path_pack = pwd()
    path_lib = "/Users/wojciechstark/Desktop/H2_on_Cu/1_h2cu_pes"
    #ush!(LOAD_PATH, path_lib)
    

    cd(path_lib) do
        
        ccall((:init_pes_, "nn.dylib"),Cvoid ,())
        
        # eann_out variables
        atom_dim = 3 # number of coordinates
        energy_out = []
        force_out = []
        atom_num = 0
        force_num = 0

        # get atoms number
        open(input_path) do file
            ln = readline(file)
            atom_num = parse(Int, ln)
        end

        force_num=atom_num-constraint_num # number of forces to calculate (number of unconstrained atoms)
        fcoor=zeros(Float64, atom_dim, atom_num) # coordinates array
        force=zeros(Float64, atom_dim, force_num) # forces array
        energy=zeros(Float64, 1, 1) # energy value (in 1x1 array because this specific fortran function requires an array not a value)
        


        open(input_path) do file
            line_num = 0 # current line overall
            subline_num = 0 # current line for current image
            points_num = 0 # current number of points calculated

            # the algorithm below is adjusted to .xyz files
            for ln in eachline(file)
                calc_pes = false
                subline_num += 1
                if subline_num > 2
                    line_text = ln
                    m = split(line_text)
                    # import a line with coordinates from the input file 
                    for k in 1:3
                        fcoor[k,subline_num-2]=parse(Float64, m[k+1])
                    end
                    if subline_num == atom_num + 2
                        subline_num = 0
                        calc_pes = true
                        points_num += 1
                    end
                end

                if calc_pes == true
                    print("coords:", fcoor, "\n\n")
                    # EANN PES: get energy and force from EANN to vars: energy and force
                    ccall((:eann_out_, "nn.dylib")
                            ,Cvoid 
                            ,(Ref{Int64}
                            ,Ref{Int64}
                            ,Ref{Float64}
                            ,Ptr{Float64}
                            ,Ptr{Float64} 
                            )
                            , coordinates_type, force_incl, fcoor, energy, force)

                    if force_incl == 0
                        push!(energy_out, energy[1,1])
                    else
                        force_out = vcat(reshape(force_out,size(force_out, 1),18), force)
                    end
                    print("Point no ", points_num, " loaded\n")

                    # save data point in an output file 
                    #if generate_output==true
                    #    open(output_path,"a") do io
                    #        println(io, "Point no ", points_num)
                    #        println(io, "Energy:")
                    #        println(io, energy[1,1])
                    #        println(io, "Force:")
                    #        for i = 1:size(force,2)
                    #            if i>1 println(io, "") end
                    #            for j = 1:size(force,1)
                    #                print(io, "   ", force[j,i])
                    #            end
                    #        end
                    #        println(io, "\n")
                    #    end
                    #end
                end
            end
            # EANN PES: call deallocate_all function
            ccall((:deallocate_all_, "nn.dylib"),Cvoid ,())
        end

        #push!(LOAD_PATH, path_pack)

        if force_incl == 1
            return force_out
        else
            return energy_out
        end
    end
end


#example use of function

#in_path = "H2_Cu.xyz" # must be of the same facet as the pes
in_path = "H2Cu_example.xyz"
output_path = "pes_output_h2Cu(111)"
constr_n = 0 # number of constrained atoms
coor_t = 0 # 0 for cart, 1 for frac
#f_incl = 1 # 1 for energy&force, 0 for energy only.
gen_output = false
#const mylib = "nn.dylib"

#print(coor_t)

out_f = H2_Cu_PES_out(in_path, coor_t, 1, constr_n, gen_output)

print(out_f, "\n\n")

out_e = H2_Cu_PES_out(in_path, coor_t, 0, constr_n, gen_output)

print(out_e)