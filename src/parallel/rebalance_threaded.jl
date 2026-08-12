using ChunkSplitters

mutable struct LoadBalancerNcoll
    n_cells::Int64
    n_chunks::Int64
    n_collisions::Vector{Float64}
    n_coll_total_per_chunk::Vector{Float64}
    n_coll_total::Float64
    chunked_indices::Vector{UnitRange{Int64}}

    function LoadBalancerNcoll(n_cells, n_chunks)
        n_chunks > n_cells && throw(ArgumentError("n_chunks ($n_chunks) cannot exceed n_cells ($n_cells)"))
        chunked_indices = index_chunks(1:n_cells; n=n_chunks)
        return new(n_cells, n_chunks, zeros(Float64, n_cells),
                   zeros(Float64, n_chunks),
                   0.0,
                   [chunk_index for chunk_index in chunked_indices])
    end
end

@inline function update_n_collisions!(lb::LoadBalancerNcoll, chunk_id, n_collisions, cell, averaging_window)
    n_c_avg = n_collisions / averaging_window
    lb.n_collisions[cell] += n_c_avg
    lb.n_coll_total_per_chunk[chunk_id] += n_c_avg
end

function rebalance_lb!(lb::LoadBalancerNcoll)
    current_chunk = 1
    lb.n_coll_total = sum(lb.n_coll_total_per_chunk)
    colls_per_chunk = lb.n_coll_total / lb.n_chunks
    coll_counter = lb.n_collisions[1]
    start_cell = 1
    for i in 1:lb.n_cells-1
        if coll_counter + lb.n_collisions[i+1] > colls_per_chunk
            # find which one is actually closer to colls_per_chunk
            if abs(coll_counter - colls_per_chunk) < abs(coll_counter + lb.n_collisions[i+1] - colls_per_chunk)
                # current cell is end of new chunk
                coll_counter = lb.n_collisions[i+1]
                lb.chunked_indices[current_chunk] = start_cell:i
                current_chunk += 1
                start_cell = i+1
            elseif coll_counter > colls_per_chunk
                # current cell is end of new chunk
                coll_counter = lb.n_collisions[i+1]
                lb.chunked_indices[current_chunk] = start_cell:i
                current_chunk += 1
                start_cell = i+1
            else
                # postpone until next chunk
                coll_counter += lb.n_collisions[i+1]
            end
        else
            coll_counter += lb.n_collisions[i+1]
        end

        if current_chunk == lb.n_chunks
            lb.chunked_indices[current_chunk] = start_cell:lb.n_cells
            return nothing
        end

        # check whether we have space to fit remaining chunks
        if lb.n_chunks - current_chunk == lb.n_cells - i
            # need to fit remaining chunks with size 1 and exit
            lb.chunked_indices[current_chunk] = start_cell:i
            current_chunk += 1
            for j in i+1:lb.n_cells
                lb.chunked_indices[current_chunk] = j:j
                current_chunk += 1
            end
            return nothing
        end
    end
end

function reset_lb!(lb::LoadBalancerNcoll)
    fill!(lb.n_collisions, 0.0)
    fill!(lb.n_coll_total_per_chunk, 0.0)
    lb.n_coll_total = 0.0
end