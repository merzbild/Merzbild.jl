# The NNLS.jl package is licensed under the MIT "Expat" License:

# Copyright (c) 2017: Robin Deits.

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Code taken from https://github.com/rdeits/NNLS.jl
# This is a cut-down version of the NNLS.jl code as provided above, removing anything related to modules and QP
@muladd begin
using LinearAlgebra

"""
CONSTRUCTION AND/OR APPLICATION OF A SINGLE
HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B

The original version of this code was developed by
Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
1973 JUN 12, and published in the book
"SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
Revised FEB 1995 to accompany reprinting of the book by SIAM.
"""
function construct_householder!(u::AbstractVector{T}, up::T)::T where T
    m = length(u)
    if m <= 1
        return up
    end

    # four independent maxima so that the reduction is not latency-bound (the maximum is exact in any order)
    cl1 = zero(T); cl2 = zero(T); cl3 = zero(T); cl4 = zero(T)
    i = 1
    @inbounds while i + 3 <= m
        a1 = abs(u[i]); a2 = abs(u[i + 1]); a3 = abs(u[i + 2]); a4 = abs(u[i + 3])
        cl1 = ifelse(a1 > cl1, a1, cl1)
        cl2 = ifelse(a2 > cl2, a2, cl2)
        cl3 = ifelse(a3 > cl3, a3, cl3)
        cl4 = ifelse(a4 > cl4, a4, cl4)
        i += 4
    end
    @inbounds while i <= m
        a1 = abs(u[i])
        cl1 = ifelse(a1 > cl1, a1, cl1)
        i += 1
    end
    cl = max(max(cl1, cl2), max(cl3, cl4))
    # @assert cl > 0
    clinv = one(T) / cl
    sm = zero(T)
    @inbounds @simd for i in 1:m
        sm += (u[i] * clinv)^2
    end
    cl *= sqrt(sm)
    @inbounds if u[1] > zero(T)
        cl = -cl
    end
    @inbounds result = u[1] - cl
    @inbounds u[1] = cl

    return result
end

"""
CONSTRUCTION AND/OR APPLICATION OF A SINGLE
HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B

The original version of this code was developed by
Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
1973 JUN 12, and published in the book
"SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
Revised FEB 1995 to accompany reprinting of the book by SIAM.
"""
function apply_householder!(u::AbstractVector{T}, up::T, c::AbstractVector{T}) where T
    m = length(u)
    if m > 1
        # cl = abs(u1)
        # @assert cl > zero(T)
        @inbounds b = up * u[1]
        if b >= zero(T)
            return
        end
        b = one(T) / b

        @inbounds sm = c[1] * up
        @inbounds @simd for i in 2:m
            sm = sm + c[i] * u[i]
        end
        if sm != zero(T)
            sm *= b
            @inbounds c[1] = c[1] + sm * up
            @inbounds @simd for i in 2:m
                c[i] = c[i] + sm * u[i]
            end
        end
    end
end

"""
    apply_householder_sweep!(A::Matrix{T}, up::T, ju, nsetp, m, idx, iz1, iz2, b, w) where T

Apply the Householder transformation held in rows `nsetp:m` of column `ju` of `A` to
rows `nsetp:m` of every column in set Z, i.e. the columns `idx[iz1], ..., idx[iz2]`,
and compute the dual vector entries `w[idx[jz]] = A[nsetp+1:m, idx[jz]] ⋅ b[nsetp+1:m]`
of those columns from the transformed values.

This is a batched replacement for calling [`apply_householder!`](@ref) once per column
followed by the dual vector computation of the next NNLS main-loop iteration.
It is equivalent to that, but the reflector scalars (the `up * u[1]` sign test and
the reciprocal `1 / (up * u[1])`) are the same for every column, so they are computed
once outside the column loop instead of once per column, the columns are processed
four at a time so that the reflector stays in registers and the four dot products
overlap, and the dual entries are accumulated while the transformed column values are
still in registers.

# Positional arguments
* `A`: the matrix being triangularized
* `up`: the first component of the Householder vector, as returned by [`construct_householder!`](@ref)
* `ju`: index of the column holding the Householder vector
* `nsetp`: number of columns already moved to set P (the transformation acts on rows `nsetp:m`)
* `m`: number of rows of `A`
* `idx`: the index vector of the NNLS workspace
* `iz1`, `iz2`: the range of positions in `idx` holding the set-Z columns
* `b`: the (already transformed) right-hand side vector
* `w`: the dual vector

# Returns
Tuple `(computed, wmax, izmax)`: `computed` is `true` if the dual vector entries of the set-Z columns
were computed and `false` if the transformation is trivial and nothing was done; `wmax` and `izmax`
are the largest positive dual entry and its position in `idx`, as returned by [`largest_positive_dual`](@ref)
(both zero if not computed or no positive entry exists).
"""
function apply_householder_sweep!(A::Matrix{T}, up::T, ju, nsetp, m, idx, iz1, iz2, b, w) where T
    len = m - nsetp + 1
    wmax = zero(T)
    izmax = zero(iz1)
    if len <= 1
        return false, wmax, izmax
    end

    # linear index of A[nsetp, ju]; the Householder vector is A[ou], A[ou+1], ..., A[ou+len-1]
    ou = nsetp + (ju - 1) * m
    @inbounds bb = up * A[ou]
    if bb >= zero(T)
        return false, wmax, izmax
    end
    binv = one(T) / bb

    jz = iz1
    @inbounds while jz + 3 <= iz2
        o1 = nsetp + (idx[jz] - 1) * m
        o2 = nsetp + (idx[jz + 1] - 1) * m
        o3 = nsetp + (idx[jz + 2] - 1) * m
        o4 = nsetp + (idx[jz + 3] - 1) * m

        s1 = A[o1] * up; s2 = A[o2] * up
        s3 = A[o3] * up; s4 = A[o4] * up
        @simd for i in 1:len-1
            ui = A[ou + i]
            s1 = s1 + A[o1 + i] * ui
            s2 = s2 + A[o2 + i] * ui
            s3 = s3 + A[o3 + i] * ui
            s4 = s4 + A[o4 + i] * ui
        end
        s1 *= binv; s2 *= binv
        s3 *= binv; s4 *= binv

        A[o1] = A[o1] + s1 * up
        A[o2] = A[o2] + s2 * up
        A[o3] = A[o3] + s3 * up
        A[o4] = A[o4] + s4 * up
        w1 = zero(T); w2 = zero(T)
        w3 = zero(T); w4 = zero(T)
        @simd ivdep for i in 1:len-1
            ui = A[ou + i]
            bi = b[nsetp + i]
            a1 = A[o1 + i] + s1 * ui
            a2 = A[o2 + i] + s2 * ui
            a3 = A[o3 + i] + s3 * ui
            a4 = A[o4 + i] + s4 * ui
            A[o1 + i] = a1
            A[o2 + i] = a2
            A[o3 + i] = a3
            A[o4 + i] = a4
            w1 = w1 + a1 * bi
            w2 = w2 + a2 * bi
            w3 = w3 + a3 * bi
            w4 = w4 + a4 * bi
        end
        w[idx[jz]] = w1
        w[idx[jz + 1]] = w2
        w[idx[jz + 2]] = w3
        w[idx[jz + 3]] = w4
        if w1 > wmax
            wmax = w1; izmax = jz
        end
        if w2 > wmax
            wmax = w2; izmax = jz + 1
        end
        if w3 > wmax
            wmax = w3; izmax = jz + 2
        end
        if w4 > wmax
            wmax = w4; izmax = jz + 3
        end
        jz += 4
    end

    # remainder, fewer than 4 columns left
    @inbounds while jz <= iz2
        oc = nsetp + (idx[jz] - 1) * m
        sm = A[oc] * up
        @simd for i in 1:len-1
            sm = sm + A[oc + i] * A[ou + i]
        end
        sm *= binv
        A[oc] = A[oc] + sm * up
        wc = zero(T)
        @simd ivdep for i in 1:len-1
            ac = A[oc + i] + sm * A[ou + i]
            A[oc + i] = ac
            wc = wc + ac * b[nsetp + i]
        end
        w[idx[jz]] = wc
        if wc > wmax
            wmax = wc; izmax = jz
        end
        jz += 1
    end
    return true, wmax, izmax
end

"""
   COMPUTE ORTHOGONAL ROTATION MATRIX..
The original version of this code was developed by
Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
1973 JUN 12, and published in the book
"SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
Revised FEB 1995 to accompany reprinting of the book by SIAM.

   COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
                      (-S,C)         (-S,C)(B)   (   0          )
   COMPUTE SIG = SQRT(A**2+B**2)
      SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
      SIG MAY BE IN THE SAME LOCATION AS A OR B .
"""
function orthogonal_rotmat(a::T, b::T)::Tuple{T, T, T} where T
    if abs(a) > abs(b)
        xr = b / a
        yr = sqrt(1.0 + xr^2)
        c = (1.0 / yr) * copysign(one(T), a)
        s = c * xr
        sig = abs(a) * yr
    elseif b != zero(T)
        xr = a / b
        yr = sqrt(1.0 + xr^2)
        s = (1.0 / yr) * copysign(one(T), b)
        c = s * xr
        sig = abs(b) * yr
    else
        sig = zero(T)
        c = zero(T)
        s = one(T)
    end
    return c, s, sig
end

"""
The original version of this code was developed by
Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
1973 JUN 15, and published in the book
"SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
Revised FEB 1995 to accompany reprinting of the book by SIAM.

Column-oriented back substitution processing four columns per pass over `zz`, every entry of `zz`
receives the same sequence of operations as in the one-column-at-a-time substitution.
"""
function solve_triangular_system!(zz, A, idx, nsetp, jj)
    q = nsetp
    @inbounds while q >= 4
        c1 = idx[q]
        c2 = idx[q - 1]
        c3 = idx[q - 2]
        c4 = idx[q - 3]

        z1 = zz[q] / A[q, c1]
        t2 = zz[q - 1] - A[q - 1, c1] * z1
        z2 = t2 / A[q - 1, c2]
        t3 = zz[q - 2] - A[q - 2, c1] * z1
        t3 = t3 - A[q - 2, c2] * z2
        z3 = t3 / A[q - 2, c3]
        t4 = zz[q - 3] - A[q - 3, c1] * z1
        t4 = t4 - A[q - 3, c2] * z2
        t4 = t4 - A[q - 3, c3] * z3
        z4 = t4 / A[q - 3, c4]

        zz[q] = z1
        zz[q - 1] = z2
        zz[q - 2] = z3
        zz[q - 3] = z4

        @simd ivdep for ii in 1:q-4
            t = zz[ii] - A[ii, c1] * z1
            t = t - A[ii, c2] * z2
            t = t - A[ii, c3] * z3
            zz[ii] = t - A[ii, c4] * z4
        end
        q -= 4
    end

    @inbounds while q >= 1
        jj = idx[q]
        zm = zz[q] / A[q, jj]
        zz[q] = zm
        for ii in 1:q-1
            zz[ii] = zz[ii] - A[ii, jj] * zm
        end
        q -= 1
    end

    if nsetp >= 1
        @inbounds jj = idx[1]
    end
    return jj
end

mutable struct NNLSWorkspace{T, I <: Integer}
    QA::Matrix{T}
    Qb::Vector{T}
    x::Vector{T}
    w::Vector{T}
    zz::Vector{T}
    idx::Vector{I}
    rnorm::T
    mode::I
    nsetp::I

    function NNLSWorkspace{T, I}(m, n) where {T, I <: Integer}
        new{T, I}(Matrix{T}(undef, m, n), # A
            Vector{T}(undef, m),    # b
            Vector{T}(undef, n),    # x
            Vector{T}(undef, n),    # w
            Vector{T}(undef, m),    # zz
            Vector{I}(undef, n),    # idx
            zero(T), # rnorm
            zero(I), # mode
            zero(I)  # nsetp
        )
    end
end

function Base.resize!(work::NNLSWorkspace{T}, m::Integer, n::Integer) where T
    work.QA = Matrix{T}(undef, m, n)
    work.Qb = Vector{T}(undef, m)
    resize!(work.x, n)
    resize!(work.w, n)
    resize!(work.zz, m)
    resize!(work.idx, n)
end

function load!(work::NNLSWorkspace{T}, A::AbstractMatrix{T}, b::AbstractVector{T}) where T
    m, n = size(A)
    @assert size(b) == (m,)
    if size(work.QA, 1) != m || size(work.QA, 2) != n
        resize!(work, m, n)
    end
    work.QA .= A
    work.Qb .= b
    work
end

NNLSWorkspace(m::Integer, n::Integer,
              eltype::Type{T}=Float64,
              indextype::Type{I}=Int) where {T, I} = NNLSWorkspace{T, I}(m, n)

function NNLSWorkspace(A::Matrix{T}, b::Vector{T}, indextype::Type{I}=Int) where {T, I}
    m, n = size(A)
    @assert size(b) == (m,)
    work = NNLSWorkspace{T, I}(m, n)
    load!(work, A, b)
    work
end


"""
Views in Julia still allocate some memory (since they need to keep
a reference to the original array). This type allocates no memory
and does no bounds checking. Use it with caution.
"""
struct UnsafeVectorView{T} <: AbstractVector{T}
    offset::Int
    len::Int
    ptr::Ptr{T}
end

UnsafeVectorView(parent::DenseArray{T}, start_ind::Integer, len::Integer) where {T} = UnsafeVectorView{T}(start_ind - 1, len, pointer(parent))
Base.size(v::UnsafeVectorView) = (v.len,)
Base.getindex(v::UnsafeVectorView, idx) = unsafe_load(v.ptr, idx + v.offset)
Base.setindex!(v::UnsafeVectorView, value, idx) = unsafe_store!(v.ptr, value, idx + v.offset)
Base.length(v::UnsafeVectorView) = v.len
Base.IndexStyle(::Type{V}) where {V <: UnsafeVectorView} = Base.IndexLinear()

"""
UnsafeVectorView only works for isbitstype types. For other types, we're already
allocating lots of memory elsewhere, so creating a new View is fine.

This function looks type-unstable, but the isbitstype(T) test can be evaluated
by the compiler, so the result is actually type-stable.
"""
@inline function fastview(parent::Array{T}, start_ind::Integer, len::Integer) where T
    if isbitstype(T)
        UnsafeVectorView(parent, start_ind, len)
    else
        @view(parent[start_ind:(start_ind + len - 1)])
    end
end

"""
Fallback for non-contiguous arrays, for which UnsafeVectorView does not make
sense.
"""
@inline fastview(parent::AbstractArray, start_ind::Integer, len::Integer) = @view(parent[start_ind:(start_ind + len - 1)])

@noinline function checkargs(work::NNLSWorkspace)
    m, n = size(work.QA)
    @assert size(work.Qb) == (m,)
    @assert size(work.x) == (n,)
    @assert size(work.w) == (n,)
    @assert size(work.zz) == (m,)
    @assert size(work.idx) == (n,)
end

function largest_positive_dual(w::AbstractVector{T}, idx::AbstractVector{TI}, s, e) where {T, TI}
    wmax = zero(T)
    izmax = zero(TI)
    @inbounds for i in s:e
        j = idx[i]
        if w[j] > wmax
            wmax = w[j]
            izmax = i
        end
    end
    wmax, izmax
end


"""
    solve!(work::NNLSWorkspace{T, TI}, ncols::Integer, max_iter::Integer, kkt_tol) where {T, TI}

Algorithm NNLS: NONNEGATIVE LEAST SQUARES

The original version of this code was developed by
Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
1973 JUN 15, and published in the book
"SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
Revised FEB 1995 to accompany reprinting of the book by SIAM.

GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM
                 A * X = B  SUBJECT TO X .GE. 0

Only the first `ncols` columns of `work.QA` are used, the workspace may hold more columns.
"""
function solve!(work::NNLSWorkspace{T, TI}, ncols::Integer, max_iter::Integer, kkt_tol) where {T, TI}
    # checkargs(work)

    A = work.QA
    b = work.Qb
    x = work.x
    w = work.w
    zz = work.zz
    idx = work.idx
    factor = 0.01
    work.mode = 1

    m = convert(TI, size(A, 1))
    n = convert(TI, ncols)

    iter = 0
    
    for i in Base.OneTo(n)
        idx[i] = i
    end

    fill!(x, zero(T))
    fill!(w, zero(T))

    # idx .= 1:n

    iz2 = n
    iz1 = one(TI)
    iz = zero(TI)
    j = zero(TI)
    jj = zero(TI)
    nsetp = zero(TI)
    up = zero(T)

    terminated = false
    # the Householder sweep computes the dual vector entries of the set-Z columns
    # and their maximum from the transformed values, so they need not be recomputed
    # unless the secondary loop modified A and b
    w_valid = false
    wmax = zero(T)
    izmax = zero(TI)

    # ******  MAIN LOOP BEGINS HERE  ******
    while true
        # println("jl main loop")
        # QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
        # OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.
        if (iz1 > iz2 || nsetp >= m)
            terminated = true
            break
        end

        # COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
        if !w_valid
            @inbounds for i in iz1:iz2
                idxi = idx[i]
                sm = zero(T)
                @simd for l in (nsetp + 1):m
                    sm = sm + A[l, idxi] * b[l]
                end
                w[idxi] = sm
            end
        end

        @inbounds while true
            # FIND LARGEST POSITIVE W(J).
            if !w_valid
                wmax, izmax = largest_positive_dual(w, idx, iz1, iz2)
            end
            w_valid = false

            # IF WMAX .LE. 0. GO TO TERMINATION.
            # THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
            if wmax <= kkt_tol
                terminated = true
                break
            end

            iz = izmax
            @inbounds j = idx[iz]

            # THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
            # BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
            # NEAR LINEAR DEPENDENCE.
            Asave = A[nsetp + 1, j]
            unorm::T = zero(T)
            @inbounds for l in 1:nsetp
                unorm += A[l, j]^2
            end
            unorm = sqrt(unorm)
            up = construct_householder!(
                # fastview(A, LinearIndices(A)[nsetp + 1, j], m - nsetp),
                fastview(A, nsetp + 1 + (j-1)*m, m - nsetp),
                up)

            @inbounds if ((unorm + abs(A[nsetp + 1, j]) * factor) - unorm) > zero(T)
                # COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
                # AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).
                # println("copying b into zz")
                # @inbounds for iii in Base.OneTo(m)
                #     zz[iii] = b[iii]
                # end
                unsafe_copyto!(zz, 1, b, 1, m)
                # zz .= b
                
                apply_householder!(
                    fastview(A, nsetp + 1 + (j-1)*m, m - nsetp),
                    up,
                    fastview(zz, nsetp + 1, m - nsetp))
                @inbounds ztest = zz[nsetp + 1] / A[nsetp + 1, j]

                # SEE IF ZTEST IS POSITIVE
                if ztest > zero(T)
                    break
                end
            end

            # REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
            # RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
            # COEFFS AGAIN.
            @inbounds A[nsetp + 1, j] = Asave
            @inbounds w[j] = zero(T)
        end
        if terminated
            break
        end

        # THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
        # SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
        # TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
        # COL J,  SET W(J)=0.
        # @inbounds for iii in Base.OneTo(m)
        #     b[iii] = zz[iii]
        # end
        unsafe_copyto!(b, 1, zz, 1, m)

        @inbounds idx[iz] = idx[iz1]
        @inbounds idx[iz1] = j
        iz1 += one(TI)
        nsetp += one(TI)

        if iz1 <= iz2
            w_valid, wmax, izmax = apply_householder_sweep!(A, up, j, nsetp, m, idx, iz1, iz2, b, w)
        end

        if nsetp != m
            @inbounds for l in (nsetp + 1):m
                A[l, j] = zero(T)
            end
        end

        @inbounds w[j] = zero(T)

        # SOLVE THE TRIANGULAR SYSTEM.
        # STORE THE SOLUTION TEMPORARILY IN ZZ().
        jj = solve_triangular_system!(zz, A, idx, nsetp, jj)

        # ******  SECONDARY LOOP BEGINS HERE ******
        #
        # ITERATION COUNTER.
        @inbounds while true
            iter += 1
            if iter > max_iter
                work.mode = 3
                terminated = true
                # println("NNLS quitting on iteration count")
                break
            end

            # SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.
            # IF NOT COMPUTE ALPHA.
            alpha = convert(T, 2)
            @inbounds for ip in Base.OneTo(nsetp)
                l = idx[ip]
                if zz[ip] <= zero(T)
                    t = -x[l] / (zz[ip] - x[l])
                    if alpha > t
                        alpha = t
                        jj = ip
                    end
                end
            end

            # IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL
            # STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.
            if alpha == 2
                break
            end

            # A AND B ARE ABOUT TO BE MODIFIED, THE DUAL VECTOR MUST BE RECOMPUTED
            w_valid = false

            # OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0 AND 1 TO
            # INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.
            @inbounds for ip in Base.OneTo(nsetp)
                l = idx[ip]
                x[l] = x[l] + alpha * (zz[ip] - x[l])
            end

            # MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I
            # FROM SET P TO SET Z.
            @inbounds i = idx[jj]

            while true
                @inbounds x[i] = zero(T)

                if jj != nsetp
                    jj += one(TI)
                    @inbounds ii = idx[jj]
                    @inbounds cc, ss, sig = orthogonal_rotmat(A[jj - 1, ii], A[jj, ii])
                    @inbounds for j in jj:nsetp
                        ii = idx[j]
                        idx[j - 1] = ii
                        A[j - 1, ii] = sig
                        A[j, ii] = zero(T)

                        # Apply procedure G2 (CC,SS,A(J-1,L),A(J,L)) to the columns other than II. The set-P columns at
                        # positions below J - 1 are zero in rows J - 1 and J and are skipped. The column at position J + 1
                        # is rotated first and the next rotation computed from it, to overlap its latency with the rest
                        p_start = j + 1
                        cc_next = cc
                        ss_next = ss
                        sig_next = sig
                        if j < nsetp
                            l = idx[j + 1]
                            temp = A[j - 1, l]
                            t2 = A[j, l]
                            A[j - 1, l] = cc * temp + ss * t2
                            A[j, l] = -ss * temp + cc * t2
                            cc_next, ss_next, sig_next = orthogonal_rotmat(A[j, l], A[j + 1, l])
                            p_start = j + 2
                        end
                        for p in p_start:n
                            l = idx[p]
                            temp = A[j - 1, l]
                            t2 = A[j, l]
                            A[j - 1, l] = cc * temp + ss * t2
                            A[j, l] = -ss * temp + cc * t2
                        end
                        temp = A[j - 1, i]
                        t2 = A[j, i]
                        A[j - 1, i] = cc * temp + ss * t2
                        A[j, i] = -ss * temp + cc * t2

                        # Apply procedure G2 (CC,SS,B(J-1),B(J))
                        temp = b[j - 1]
                        b[j - 1] = cc * temp + ss * b[j]
                        b[j] = -ss * temp + cc * b[j]

                        cc = cc_next
                        ss = ss_next
                        sig = sig_next
                    end
                end

                nsetp -= one(TI)
                iz1 -= one(TI)
                @inbounds idx[iz1] = i

                # SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
                # BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
                # IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
                # THAT ARE NONPOSITIVE WILL BE SET TO ZERO
                # AND MOVED FROM SET P TO SET Z.
                allfeasible = true
                @inbounds for jji in Base.OneTo(nsetp)
                    i = idx[jji]
                    if x[i] <= zero(T)
                        allfeasible = false
                        jj = jji
                        break
                    end
                end
                if allfeasible
                    break
                end
            end

            # COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
            # zz .= b
            # @inbounds for iii in Base.OneTo(m)
            #     zz[iii] = b[iii]
            # end
            unsafe_copyto!(zz, 1, b, 1, m)
            jj = solve_triangular_system!(zz, A, idx, nsetp, jj)
        end
        if terminated
            break
        end
        # ******  END OF SECONDARY LOOP  ******

        @inbounds for i in 1:nsetp
            x[idx[i]] = zz[i]
        end
        # ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.
    end

    # ******  END OF MAIN LOOP  ******
    # COME TO HERE FOR TERMINATION.
    # COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.
    # println("ITERATION COUNT $iter")
    sm = zero(T)
    if nsetp < m
        for i in (nsetp + 1):m
            sm += b[i]^2
        end
    # else
    #     w .= 0
    end
    work.rnorm = sqrt(sm)
    work.nsetp = nsetp
    return work.x
end

"""
    solve!(work::NNLSWorkspace, max_iter::Integer=(3 * size(work.QA, 2)), kkt_tol=1e-16)

Solve the NNLS problem posed by all columns of `work.QA`.
"""
solve!(work::NNLSWorkspace, max_iter::Integer=(3 * size(work.QA, 2)), kkt_tol=1e-16) = solve!(work, size(work.QA, 2), max_iter, kkt_tol)

function solve!(work::NNLSWorkspace{T}, A::AbstractMatrix{T}, b::AbstractVector{T}, max_iter=(3 * size(A, 2))) where T
    load!(work, A, b)
    solve!(work, max_iter)
    work.x
end

function nnls(A::DenseMatrix{T}, b::DenseVector{T}, max_iter=(3 * size(A, 2))) where T
    work = NNLSWorkspace(A, b)
    solve!(work, max_iter)
    work.x
end
end
