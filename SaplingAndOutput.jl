function Z_iteration(L)
    num_iteration = 10e4
    partition_array = []
    for i = 1:num_iteration
        partition_per_sample = partition_recursive_calculation(L)
        push!(partition_array, partition_per_sample)
    end
    mean_value = mean(partition_array)
    error_bar = std(partition_array) / sqrt(num_iteration)
    return mean_value, error_bar
end

function write_to_file(mean_value, error_bar, Δτ)
    output_file = open("result","a+")
    write(output_file, string(β), " ", string(Δτ), " ", string(mean_value), " ", string(error_bar), "\n")
    close(output_file)
end
