
function add_text_blocks(file)
    preamble = """
    ```@setup logging
    @info "Expanding $file..."
    start_time = time()
    ```
    """

    postamble = """
    ```@setup logging
    runtime = round(time() - start_time; digits=2)
    @info "...done after \$runtime s."
    ```
    """

    (tmppath, tmpio) = mktemp()
    open(file) do io
        write(tmpio, preamble)
        for line in eachline(io, keep=true) # keep so the new line isn't chomped
            write(tmpio, line)
        end
        write(tmpio, postamble)
    end
    close(tmpio)
    mv(tmppath, file, force=true)
end


for (root, dirs, files) in walkdir("src")
    for file in files
        if splitext(file)[2] == ".md"
            filepath = joinpath(root, file)
            add_text_blocks(filepath)
        end
    end
end
