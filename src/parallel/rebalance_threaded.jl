using ChunkSplitters

"""
    LoadBalancerCellQ

Used to perform load-balancing in multi-threaded simulations by tracking a per-cell
quantity (i.e. number of collisions, number of particles, density, etc.) in a range of cells indexed
by `cell_start:cell_end`
and allowing to rebalance the assignment of cells to chunks based on the total quantity per chunk.

# Fields
* `cell_start`: index of first cell
* `cell_end`: index of last cell
* `n_cells`: number of cells
* `n_chunks`: number of chunks
* `q`: `Vector` of length `n_cells` of per-cell quantity on which balancing is based. Note that this array is 1-indexed.
* `q_total_per_chunk`: `Vector` of length `n_chunks` holding the sums of `q` across all cells in a chunk
* `q_total`: total sum of the per-cell quantity across all `n_cells` cells
* `chunked_indices`: a `Vector{UnitRange{Int64}}` of length `n_chunks` holding the ranges of indices
assigned to each chunk
"""
mutable struct LoadBalancerCellQ
    cell_start::Int64
    cell_end::Int64
    n_cells::Int64
    n_chunks::Int64
    q::Vector{Float64}
    q_total_per_chunk::Vector{Float64}
    q_total::Float64
    chunked_indices::Vector{UnitRange{Int64}}

    @doc """
        LoadBalancerCellQ(n_cells, n_chunks)

    Create a per-cell indicator-based load balancer for `n_cells` cells and `n_chunks` chunks.
    Note: this uses `index_chunks`, i.e. it is assumed that cell indexing starts from 1 and is contiguous!

    # Positional arguments
    * `n_cells`: number of cells
    * `n_chunks`: number of chunks

    # Keyword arguments
    * `cell_start`: index of first cell indexed by the `LoadBalancerCellQ` instance
    """
    function LoadBalancerCellQ(n_cells, n_chunks; cell_start=1)
        n_chunks > n_cells && throw(ArgumentError("n_chunks ($n_chunks) cannot exceed n_cells ($n_cells)"))
        cell_end = cell_start + n_cells - 1
        chunked_indices = index_chunks(cell_start:cell_end; n=n_chunks)
        return new(cell_start, cell_end, n_cells, n_chunks, zeros(Float64, n_cells),
                   zeros(Float64, n_chunks),
                   0.0,
                   [chunk_index for chunk_index in chunked_indices])
    end
end

"""
    update_lb_cellq!(lb::LoadBalancerCellQ, chunk_id, cell, q; averaging_window=1.0)

Update the values of the load-balancing quantity being tracked in a single cell.

# Positional arguments
* `lb`: the `LoadBalancerCellQ` instance to update
* `chunk_id`: the index of the chunk to which the cell belongs
* `cell`: the index of the cell to update
* `q`: the value of the load-balancing quantity in the cell

# Keyword arguments
* `averaging_window`: the time window over which the load-balancing quantity is averaged
(to potentially avoid round-off errors, overflow, or variable time-steps)
"""
@inline function update_lb_cellq!(lb::LoadBalancerCellQ, chunk_id, cell, q; averaging_window=1.0)
    n_c_avg = q / averaging_window
    lb.q[cell - lb.cell_start + 1] += n_c_avg
    lb.q_total_per_chunk[chunk_id] += n_c_avg
end

"""
    rebalance_lb!(lb::LoadBalancerCellQ)

Compute new ranges of cell-indices so that each range holds approximately the same total
sum of the load-balancing quantity.

# Positional arguments
* `lb`: the `LoadBalancerCellQ` instance on which to perform re-balancing
"""
function rebalance_lb!(lb::LoadBalancerCellQ)
    current_chunk = 1
    offset = lb.cell_start - 1
    lb.q_total = sum(lb.q_total_per_chunk)
    q_per_chunk = lb.q_total / lb.n_chunks
    q_counter = lb.q[lb.cell_start-offset]
    start_cell = lb.cell_start
    for i in lb.cell_start:lb.cell_end-1
        if q_counter + lb.q[i-offset+1] > q_per_chunk
            # find which one is actually closer to q_per_chunk
            if abs(q_counter - q_per_chunk) < abs(q_counter + lb.q[i-offset+1] - q_per_chunk)
                # current cell is end of new chunk
                q_counter = lb.q[i-offset+1]
                lb.chunked_indices[current_chunk] = start_cell:i
                current_chunk += 1
                start_cell = i+1
            elseif q_counter > q_per_chunk
                # current cell is end of new chunk
                q_counter = lb.q[i-offset+1]
                lb.chunked_indices[current_chunk] = start_cell:i
                current_chunk += 1
                start_cell = i+1
            else
                # postpone until next chunk
                q_counter += lb.q[i-offset+1]
            end
        else
            q_counter += lb.q[i-offset+1]
        end

        if current_chunk == lb.n_chunks
            lb.chunked_indices[current_chunk] = start_cell:lb.cell_end
            return nothing
        end

        # check whether we have space to fit remaining chunks
        if lb.n_chunks - current_chunk == lb.cell_end - i
            # need to fit remaining chunks with size 1 and exit
            lb.chunked_indices[current_chunk] = start_cell:i
            current_chunk += 1
            for j in i+1:lb.cell_end
                lb.chunked_indices[current_chunk] = j:j
                current_chunk += 1
            end
            return nothing
        end
    end
end

"""
    reset_lb!(lb::LoadBalancerCellQ)

Set all values of the load balancer tracking a per-cell quantity to zero.
Indexing is not reset.

# Positional arguments
* `lb`: the `LoadBalancerCellQ` instance to reset
"""
function reset_lb!(lb::LoadBalancerCellQ)
    fill!(lb.q, 0.0)
    fill!(lb.q_total_per_chunk, 0.0)
    lb.q_total = 0.0
end