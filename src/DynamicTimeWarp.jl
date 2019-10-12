module DynamicTimeWarp


"""
    dtw(
        sequence1, sequence2,
        dist=(a, b)->(a - b)'*(a - b)
    ) -> (path1, path2, cost)

Compute the set of indices `path1, path2` that give the alignment of two
sequences over their last dimension by dynamic time warping.

Accepts as an optional parameter the distance function to be used in the
evaluation of the optimum path. The two arguments passed to the distance
function are views of the original sequences over the last dimension preserving
their shape. In a vector, for example, the distances will be computed between
two scalars, in a 2 by n array between 2 by 1 vectors, in a 2 by 3 by n between
2 by 3 arrays and so on.

Returns the warping paths and the total cost of the optimum path.

See also: [`linkseq`](@ref)
"""
function dtw(
	sequence1::AbstractArray,
	sequence2::AbstractArray,
	dist::Function=(a, b)->(a - b)'*(a - b)
)

	length1 = size(sequence1)[end]
	length2 = size(sequence2)[end]


	# The DTW is calculated over the last dimension of the array, the next
	# function ensures that the shape of the original elements in the sequence
	# will be preserved during the computation of the distance.
	d = (a, b) -> dist(
		selectdim(sequence1, ndims(sequence1), a),
		selectdim(sequence2, ndims(sequence2), b)
	)

	initialCost = d(1, 1)
	cost = Matrix{typeof(initialCost)}(undef, length1, length2)
	cost[1, 1] = initialCost


	# Fill the first line...
	for i in 2:length1
		cost[i, 1] = cost[i - 1, 1] + d(i, 1)
	end

	# ...and the first column of the cost matrix.
	for j in 2:length2
		cost[1, j] = cost[1, j - 1] + d(1, j)
	end

	# Compute the accumulated cost of the remaining elements
	for i in 2:length1
		for j in 2:length2
			cost[i, j] = d(i, j) +
				min(cost[i - 1, j], cost[i - 1, j - 1], cost[i, j - 1])
		end
	end


	return trackback(cost)..., cost[end]

end


"""
    trackback(cost::AbstractMatrix) -> (path1, path2)

Compute the lowest cost path from a cost array. Returns two vectors
representing the path.
"""
function trackback(cost::AbstractMatrix)

	i, j = size(cost)
	# Start to create the path from the last element of the cost matrix.
	path1, path2 = [i], [j]

	while i > 1 && j > 1
		least = argmin(cost[i-1:i,j-1:j])
		# Everytime that the minimal cost is reached by going back along one
		# dimension, decrease the corresponding index.
		push!(path1, (least[1] == 1 ? i-=1 : i))
		push!(path2, (least[2] == 1 ? j-=1 : j))
	end

	# One of the indexes may still be greater than one if the path hits a border
	# of the cost matrix, in this case the path must be completed with a straight
	# line to the first element of the cost matrix.
	if i > 1
		append!(path1, i-1:-1:1)
		append!(path2, ones(length(path1) - length(path2)))
	end

	if j > 1
		append!(path2, j-1:-1:1)
		append!(path1, ones(length(path2) - length(path1)))
	end


	return reverse(path1), reverse(path2)
end


"""
    linkseq(sequence1, sequence2, sr=1, noflinks=0) -> (x, y)

Create vectors containing the horizontal and the vertical coordinates of each
pair of corresponding points of two sequences computated by using the dtw.

The returned vectors can be used to plot the paths linking the two sequences
over an existing plot:
```jldoctest
julia> p = plot(0:1/100:length(sequence1), [sequence1, sequence2])
julia> x, y = linkseq(sequence1, sequence2, 100, 20)
julia> plot!(p, x, y);
```


**Parameters:**

`sequence1`, `sequence2`: Sequences for evaluating the dtw

`sr`: The sampling rate to be used in the horizontal axis, should match the
original sampling rate of the sequence.

`noflinks`: The number of links between pairs of corresponding points to be
created. If zero is received create links between all pairs.
"""
function linkseq(sequence1::AbstractArray, sequence2::AbstractArray, sr=1, noflinks=0)

	path1, path2 = dtw(sequence1, sequence2)

	# Recreate the original independent axis
	len1 = size(sequence1)[end]
	len2 = size(sequence2)[end]
	t = range(0; step=1/sr, length=max(len1, len2))

	links = 1 : (noflinks == 0 ? 1 : (length(path1) ÷ noflinks)) : length(path1)

	# Create an array of ranges that will be used to preserve the original shape
	# of the sequence
	s = map(n->1:n, size(sequence1)[1:end-1])

	x = Array{typeof(t[1])}(undef, 2, length(links))
	y = Array{typeof(sequence1[1])}(undef, 2, length.(s)..., length(links))

	j = 1
	for l in links
		x[:,j] = [t[path1[l]]; t[path2[l]]]
		y[:, s..., j] = [sequence1[s..., path1[l]]; sequence2[s..., path2[l]]]
		j+=1
	end

	return (x, y)

end


export dtw, linkseq

end # module
