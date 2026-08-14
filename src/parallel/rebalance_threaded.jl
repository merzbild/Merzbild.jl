using ChunkSplitters

mutable struct LoadBalancerCellQ
    n_cells::Int64
    n_chunks::Int64
    q::Vector{Float64}
    q_total_per_chunk::Vector{Float64}
    q_total::Float64
    chunked_indices::Vector{UnitRange{Int64}}

    function LoadBalancerCellQ(n_cells, n_chunks)
        n_chunks > n_cells && throw(ArgumentError("n_chunks ($n_chunks) cannot exceed n_cells ($n_cells)"))
        chunked_indices = index_chunks(1:n_cells; n=n_chunks)
        return new(n_cells, n_chunks, zeros(Float64, n_cells),
                   zeros(Float64, n_chunks),
                   0.0,
                   [chunk_index for chunk_index in chunked_indices])
    end
end

@inline function update_lb_cellq!(lb::LoadBalancerCellQ, chunk_id, q, cell, averaging_window)
    n_c_avg = q / averaging_window
    lb.q[cell] += n_c_avg
    lb.q_total_per_chunk[chunk_id] += n_c_avg
end

function rebalance_lb!(lb::LoadBalancerCellQ)
    current_chunk = 1
    lb.q_total = sum(lb.q_total_per_chunk)
    q_per_chunk = lb.q_total / lb.n_chunks
    q_counter = lb.q[1]
    start_cell = 1
    for i in 1:lb.n_cells-1
        if q_counter + lb.q[i+1] > q_per_chunk
            # find which one is actually closer to q_per_chunk
            if abs(q_counter - q_per_chunk) < abs(q_counter + lb.q[i+1] - q_per_chunk)
                # current cell is end of new chunk
                q_counter = lb.q[i+1]
                lb.chunked_indices[current_chunk] = start_cell:i
                current_chunk += 1
                start_cell = i+1
            elseif q_counter > q_per_chunk
                # current cell is end of new chunk
                q_counter = lb.q[i+1]
                lb.chunked_indices[current_chunk] = start_cell:i
                current_chunk += 1
                start_cell = i+1
            else
                # postpone until next chunk
                q_counter += lb.q[i+1]
            end
        else
            q_counter += lb.q[i+1]
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

function reset_lb!(lb::LoadBalancerCellQ)
    fill!(lb.q, 0.0)
    fill!(lb.q_total_per_chunk, 0.0)
    lb.q_total = 0.0
end