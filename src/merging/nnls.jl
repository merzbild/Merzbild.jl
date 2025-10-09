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

    cl = maximum(abs, u)
    @assert cl > 0
    clinv = 1 / cl
    sm = zero(T)
    for ui in u
        sm += (ui * clinv)^2
    end
    cl *= sqrt(sm)
    @inbounds if u[1] > 0
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
        @inbounds cl = abs(u[1])
        @assert cl > 0
        @inbounds b = up * u[1]
        if b >= 0
            return
        end
        b = 1 / b

        @inbounds sm = c[1] * up
        @inbounds for i in 2:m
            sm = sm + c[i] * u[i]
        end
        if sm != 0
            sm *= b
            @inbounds c[1] = c[1] + sm * up
            @inbounds for i in 2:m
                c[i] = c[i] + sm * u[i]
            end
        end
    end
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
        yr = sqrt(1 + xr^2)
        c = (1 / yr) * sign(a)
        s = c * xr
        sig = abs(a) * yr
    elseif b != 0
        xr = a / b
        yr = sqrt(1 + xr^2)
        s = (1 / yr) * sign(b)
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
"""
function solve_triangular_system!(zz, A, idx, nsetp, jj)
    @inbounds for l in Base.OneTo(nsetp)
        ip = nsetp + 1 - l
        if (l != 1)
            for ii in 1:ip
                zz[ii] = zz[ii] - A[ii, jj] * zz[ip + 1]
            end
        end
        jj = idx[ip]
        zz[ip] /= A[ip, jj]
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
Algorithm NNLS: NONNEGATIVE LEAST SQUARES

The original version of this code was developed by
Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
1973 JUN 15, and published in the book
"SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
Revised FEB 1995 to accompany reprinting of the book by SIAM.

GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM
                 A * X = B  SUBJECT TO X .GE. 0
"""
function solve!(work::NNLSWorkspace{T, TI}, max_iter::Integer=(3 * size(work.QA, 2)), kkt_tol=1e-16) where {T, TI}
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
    n = convert(TI, size(A, 2))

    iter = 0
    
    for i in Base.OneTo(n)
        x[i] = 0.0
        idx[i] = i
        w[i] = 0.0
    end

    # idx .= 1:n

    iz2 = n
    iz1 = one(TI)
    iz = zero(TI)
    j = zero(TI)
    jj = zero(TI)
    nsetp = zero(TI)
    up = zero(T)

    terminated = false

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
        @inbounds for i in iz1:iz2
            idxi = idx[i]
            sm = zero(T)
            for l in (nsetp + 1):m
                sm = sm + A[l, idxi] * b[l]
            end
            w[idxi] = sm
        end

        @inbounds while true
            # FIND LARGEST POSITIVE W(J).
            wmax, izmax = largest_positive_dual(w, idx, iz1, iz2)

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
            up = construct_householder!(
                # fastview(A, LinearIndices(A)[nsetp + 1, j], m - nsetp),
                fastview(A, nsetp + 1 + (j-1)*m, m - nsetp),
                up)
            unorm::T = zero(T)
            @inbounds for l in 1:nsetp
                unorm += A[l, j]^2
            end
            unorm = sqrt(unorm)

            @inbounds if ((unorm + abs(A[nsetp + 1, j]) * factor) - unorm) > 0
                # COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
                # AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).
                # println("copying b into zz")
                @inbounds for iii in Base.OneTo(m)
                    zz[iii] = b[iii]
                end
                # zz .= b
                
                apply_householder!(
                    fastview(A, nsetp + 1 + (j-1)*m, m - nsetp),
                    up,
                    fastview(zz, nsetp + 1, m - nsetp))
                @inbounds ztest = zz[nsetp + 1] / A[nsetp + 1, j]

                # SEE IF ZTEST IS POSITIVE
                if ztest > 0
                    break
                end
            end

            # REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
            # RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
            # COEFFS AGAIN.
            @inbounds A[nsetp + 1, j] = Asave
            @inbounds w[j] = 0
        end
        if terminated
            break
        end

        # THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
        # SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
        # TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
        # COL J,  SET W(J)=0.
        @inbounds for iii in Base.OneTo(m)
            b[iii] = zz[iii]
        end

        @inbounds idx[iz] = idx[iz1]
        @inbounds idx[iz1] = j
        iz1 += one(TI)
        nsetp += one(TI)

        if iz1 <= iz2
            @inbounds for jz in iz1:iz2
                jj = idx[jz]
                apply_householder!(
                    fastview(A, nsetp + (j-1)*m, m - nsetp + 1),
                    up,
                    fastview(A, nsetp + (jj-1)*m, m - nsetp + 1))
            end
        end

        if nsetp != m
            @inbounds for l in (nsetp + 1):m
                A[l, j] = 0
            end
        end

        @inbounds w[j] = 0

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
                if zz[ip] <= 0
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
                x[i] = 0

                if jj != nsetp
                    jj += one(TI)
                    @inbounds for j in jj:nsetp
                        ii = idx[j]
                        idx[j - 1] = ii
                        cc, ss, sig = orthogonal_rotmat(A[j - 1, ii], A[j, ii])
                        A[j - 1, ii] = sig
                        A[j, ii] = 0
                        for l in Base.OneTo(n)
                            if l != ii
                                # Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))
                                temp = A[j - 1, l]
                                A[j - 1, l] = cc * temp + ss * A[j, l]
                                A[j, l] = -ss * temp + cc * A[j, l]
                            end
                        end

                        # Apply procedure G2 (CC,SS,B(J-1),B(J))
                        temp = b[j - 1]
                        b[j - 1] = cc * temp + ss * b[j]
                        b[j] = -ss * temp + cc * b[j]
                    end
                end

                nsetp -= one(TI)
                iz1 -= one(TI)
                idx[iz1] = i

                # SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
                # BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
                # IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
                # THAT ARE NONPOSITIVE WILL BE SET TO ZERO
                # AND MOVED FROM SET P TO SET Z.
                allfeasible = true
                @inbounds for jji in Base.OneTo(nsetp)
                    # jj = jji
                    i = idx[jji]
                    if x[i] <= 0
                        allfeasible = false
                        break
                    end
                end
                if allfeasible
                    break
                end
            end

            # COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
            # zz .= b
            @inbounds for iii in Base.OneTo(m)
                zz[iii] = b[iii]
            end
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
