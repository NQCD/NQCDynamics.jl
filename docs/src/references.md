```@setup logging
@info "Expanding src/references.md..."
start_time = time()
```
# References

```@bibliography
```
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
