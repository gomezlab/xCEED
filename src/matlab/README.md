#xCEED README

These matlab scripts are provided as a fully functional proof of concept. You should be able to get a good idea of how to use them by reading `test_vCEED.m` and `test_gCEED.m`. All xCEED methods (rCEED, vCEED, and, gCEED) use distance matrices as input. Use `read_phylip_dst.m` to read in distance files having the phylip format.

As a working test, you can run `test_vCEED.m` or `test_gCEED.m` from within the src/matlab directory. 


####`rCEED.m`
	%   Computes the rmsd from indirect superimposition of two given evolutionary
	%   structures (or point patterns)
	%
	% Input:
	%    (1) distance vector of protein family 1
	%    (2) distance vector of protein family 2
	%    (3) distance vector of reference family (e.g., 16S rRNA orthologs)
	%
	% Output:
	%    (1) Normalized rmsd between protein family 1 and reference
	%    (2) Normalized rmsd between protein family 2 and reference
	%    (3) Normalized rmsd between protein family 1 and 2
	%        (which is the rmsd of indirect superimposition)
	%
	%%%

####`vCEED.m`
	%   Computes the weighted rmsd from "robust" superimposition of two given
	%   evolutionary distance sets (Huber kernel and Bi-weight kernel are supported.)
	%
	% Input:
	%    (1) distance matrix of protein family 1
	%    (2) distance matrix of protein family 2
	%
	% Output:
	%    (1) Weighted rmsd between protein family 1 and 2
	%    (2) Embedding of protein family 1
	%    (3) Embedding of protein family 2 superimposed upon 1
	%    (4) Transformation parameters
	%
	% Note:
	%    robust_procrustes.m is a function that is internally used by vCEED.m
	%
	%%%

####`gCEED.m`
	%   Aligns two point patterns embedded in Euclidean space and recovers
	%   correspondence
	% 
	% Input:
	%    dmat1obs & dmat2obs: distance matrix 1 & 2
	%    numrep: the number of iteration to perform
	%    neighbound: the size of neighborhood (used to cluster embedded 
	%        nodes/points)
	%    nodes1 & 2: the initially assigned node id's
	%    known_pairs1 & 2: node id's of known interaction pairs
	%    dbound: distance constraint between known interacting nodes, if any
	%    method: embedding method ('stochastic proximity embedment' is supported)
	% 
	% Output:
	%   (1) Coordinates of embedded evolutionary structure (protein family 1) in
	%       Euclidean space, superimposed onto that of protein family 2
	%   (2) Coordinates of embedded evolutionary structure (protein family 2)
	%   (3) node id's of predicted interacting pairs for protein family 1
	%   (4) node id's of predicted interacting pairs for protein family 2
	%
	% Note:
	%    gceedConstraint.m, gceedCostfunc.m, generate_random_rotation.m, 
	%        generate_random_translatioin.m, and gmmreg_highdim.m are functions that are 
	%        internally used by gCEED.m during the optimization process.
	%
	%%%

####`compute_closest_points.m`
	%    Assigns correspondence between two point patterns based on distances
	%    between them. In case, we have known pairs, their pairing is reserved over.
	%
	% Input:
	%    dmat: Distance matrix (d_ij is distance between x_i and y_j)
	%    known_pairs1 & 2: node id's of known pairs in point set 1 & 2
	%    
	% Output:
	%    node id's of corresponding pairs
	%
	% Note:
	%    known_pairs 1 & 2 are column vectors
	%
	%%%

####`define_neighborhood.m`
	%    Divide nodes into clusters based on the size of neighborhood
	%
	% Input:
	%    neibound: The size of neighborhood
	%    dmat: Distance between nodes
	%
	% Output:
	%    gidvec: Vector of cluster center node id
	%
	%%%

####`embed_in_space.m`
	%    Projects points into Euclidean space that best satisfies the given distance relationship
	%
	% Input:
	%    dmat: Distances between nodes
	%    method: Method for embedment (metric MDS, non-metric MDS, or SPE)
	%    sizenei: The size of neighborhood if SPE (stochastic proximity
	%        embedment method is chosen)
	%
	% Output:
	%    X: The coordinates of points in the embedded structure
	%
	% Note:
	%    spembed.m is called internally when SPE method is chosen
	%
	%%%

####`get_gauss_weight_mat.m`
	%    Converts distances into probabilities with Gaussian distribution, and
	%    normalizes each row of the distance matrix the probabilities
	%
	% Input:
	%    X & Y: Coordinates of the point pattern X & Y
	%    scale: variance for gaussian distribution
	%
	% Output:
	%    W: Normalized probability matrix (w_ij is the probability of
	%    correspondence between x_i and y_j)
	%
	% Note:
	%    get_sqdist.m is internally called to compute squared distances between
	%    points
	%
	%%%

####`get_nearnei.m`
	%    Computes the distances between two point patterns and sorts them in
	%    the increasing order
	%
	% Input:
	%    coord1 & 2: Coordinates of point pattern 1 & 2
	%
	% Output:
	%    dmat: Distances between point pattern 1 & 2
	%    dmat_sorted: 'dmat' sorted row-wise
	%    rmat: Rankings in closeness (Each row vector lists the node id's of
	%        the most close point to the farthest
	%
	%%%

####`get_vicinity_hitrate.m`
	%    Computes the percentage of points whose known interaction partner is
	%        included in its top three (or the given threshold ranking) closest
	%        points
	%
	% Input:
	%    rmat: Rankings of nodes in terms of closeness (or vicinity)
	%    bound: Threshold ranking
	%
	% Output:
	%    score: The vicinity hitrate
	%
	%%%

####`standardize_config.m`
	%    Re-centers and rescales point pattern 
	%
	% Input:
	%    X: Coordinates of points
	%
	% Output:
	%    Z: A transformed coordinates, centered at the origin and rescaled with
	%    the standard deviation
	%
	%%%

####`read_phylip_dst.m`
	%    Reads in phylip distance file
	%
	% Input: distance file name
	%
	% Output:
	%    species: a cell array of species name  
	%    dmat: a matrix of double (symmetrci)
	%
	%%%

####`draw_3d_snapshot.m`
	%    Displays point patterns and correspondences
	%
	% Input:
	%    Z1 & 2: Coordinates of point pattern 1 & 2
	%    markercolor1 & 2: Color for point pattern 1 & 2 
	%    linecolor: Color for edge drawn between corresponding pairs
	%    labels1 & 2: Text label for each point
	%    pairs1 & 2: Node id's of corresponding pairs
	%
	% Output:
	%    None, but a figure window will pop up with the 3D plot
	%
	% Note:
	%    If labels1 = [], then no label will be printed.
	%
	%%%

