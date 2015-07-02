require(network)
require(sna)


run.simulation <- function(Authors, T.it = 10, pSingle = 0.2, colType = "Clust", name = NULL) {

	nA <- length(Authors)
	t = 0
	repeat {

		t <- t + 1

		if (t == T.it) {
			break
		}

		## Select a lead author:::
		lead_A = runif(1, 1, nA)%/%1

		## Select collaborators:::
		coll <- vector("integer")
		if (runif(1, 0, 1) < pSingle) {
			coll <- c(coll, lead_A)
		} else if (colType == "Attach") {
			coll <- seek.Preferential(Authors, lead_A)
		} else if (colType == "Clust") {
			coll <- seek.Transitivity(Authors, lead_A)
		} else if (colType == "Sim") {
			coll <- seek.Similar(Authors, lead_A)
		} else {
			coll <- seek.Diff(Authors, lead_A)
		}

		# Update co-authorship:::
		Authors <- update.AxA(Authors, coll, t)

		## Create knowledge:::
		
		# check common knowledge set:
		A.List <- vector("list")
		for (i in 1:length(coll)) {
			ind <- coll[i]
			A.List <- c(A.List, list(Authors[[ind]]))
		}
		
		common.keys <- get.common.knowledge(A.List)
		 
		# determine starting posn for semantic walk:
		start_key = NULL
		if (length(common.keys) > 0) {
			c.degrees <- get.knowledge.degrees(A.List, common.keys)
			start_key <- pref.attach.choice(c.degrees)
		}
				
		# do parallel semantic walk and merge:
		
		new_knowledge <- vector("integer")
		keys <- vector("integer")
		for (i in 1:length(coll)) {
			
			A <- A.List[[i]]
	
			if (is.null(start_key)) {
				start_key <- pick_init_knowledge(A)
			}
						 
			if (is.null(start_key)) {
				next
			}
			keys <- get.a.semantic.walk(A, start_key)
		
			new_knowledge <- union(new_knowledge, keys)
		}


		## Update knowledge space:::
		Authors <- update.KxK(Authors, coll, new_knowledge, t)


		# Report:::
		
	}

	return(T)
}

pick_init_knowledge <- function(Author) {
	links <- vector("list")
	nK <- length(Author$K.Set)
	if (nK == 0) {
		return(c())
	}

	for (i in 1:nK) {
		new <- get.knowledge.ego.net(Author, Author$K.Set[i])
		if (length(new) == 0) {
			next
		}
		links <- c(links, new)
	}
	key <- pref.attach.choice(links, is.node = F)
	return(key)
}

seek.Diff <- function(Authors, lead_A, pDiff = 0) {
	colls <- vector("integer")
	colls <- c(colls, lead_A)
	A1 <- Authors[[lead_A]]
	new <- lead_A
	diff <- 0
	new_diff <- 0
	nA <- length(Authors)
	if (nA == 0) {
		return(colls)
	}

	for (i in 1:nA) {
		A2.id <- Authors[[i]]$ID
		if (A2.id == lead_A) {
			next
		}
		A2 <- Authors[[i]]
		new_diff <- compute.Kdiff(A1, A2)
		if (new_diff > diff) {
			diff <- new_diff
			new <- A2.id
		}
	}

	colls <- union(colls, new)

	return(colls)
}

seek.Similar <- function(Authors, lead_A, pSim = 0) {
	colls <- vector("integer")
	colls <- c(colls, lead_A)
	A1 <- Authors[[lead_A]]
	new <- lead_A
	sim <- 0
	new_sim <- 0
	nA <- length(Authors)
	if (nA == 0) {
		return(colls)
	}

	for (i in 1:nA) {
		A2.id <- Authors[[i]]$ID
		if (A2.id == lead_A) {
			next
		}
		A2 <- Authors[[i]]
		new_sim <- compute.Ksimilarity(A1, A2)
		if (new_sim > sim) {
			sim <- new_sim
			new <- A2.id
		}
	}

	colls <- union(colls, new)

	return(colls)
}


seek.Preferential <- function(Authors, lead_A, pAttach = 0) {
	colls <- vector("integer")
	colls <- c(colls, lead_A)
	degrees <- get.degree.distro(Authors)
	co_A <- pref.attach.choice(degrees)
	colls <- union(colls, co_A)
	if (runif(1) < pAttach) {
		co_A <- pref.attach.choice(degrees)
		colls <- union(colls, co_A)
	}
	return(colls)
}

seek.Transitivity <- function(Auothors, lead_A) {
	colls <- vector("integer")
	colls <- c(colls, lead_A)
	links <- get.connections(lead_A, Auothors)
	co_A <- lead_A
	if (length(links) == 0) {
		degrees <- get.degree.distro(Auothors)
		co_A <- pref.attach.choice(degrees)
	}
	colls <- get.a.walk(Auothors, co_A, path = colls, k = 1)
	return(colls)
}


make.author <- function(id = 0) {
	# Author Data Structure:
	
	author <- list()
	author$ID = id
	author$Pub.Count = 0

	author$A.Set = vector(mode = "integer")
	author$AxA.id = vector(mode = "integer")
	author$AxA.w = vector(mode = "integer")
	author$AxA.t = vector(mode = "integer")

	author$K.Set = vector(mode = "integer")
	author$KxK.key1 = vector(mode = "integer")
	author$KxK.key2 = vector(mode = "integer")
	author$KxK.w = vector(mode = "integer")
	author$KxK.t = vector(mode = "integer")

	return(author)
}

initialize.community <- function(AxA.data) {
	author <- make.author()
	N = max(AxA.data$to, AxA.data$from)
	Community <- vector("list", N)
	for (i in 1:N) {
		a <- author
		a$ID <- i
		Community[[i]] <- a
	}

	for (i in 1:length(AxA.data$from)) {
		from <- AxA.data$from[i]
		to <- AxA.data$to[i]
		w <- AxA.data$w[i]
		t <- AxA.data$t[i]

		if (from == 0 & to == 0) {
			next
		}
		if (from == 0) {
			Community[[to]]$Pub.Count <- Community[[to]]$Pub.Count + 1
			next
		}
		if (to == 0) {
			Community[[from]]$Pub.Count <- Community[[from]]$Pub.Count + 1
			next
		}

		Community[[from]]$Pub.Count <- Community[[from]]$Pub.Count + 1
		Community[[to]]$Pub.Count <- Community[[to]]$Pub.Count + 1

		#Add to the collaborator set the $A.List:
		Community[[from]]$A.Set <- union(Community[[from]]$A.Set, to)
		Community[[to]]$A.Set <- union(Community[[to]]$A.Set, from)


		Community[[from]]$AxA.id <- c(Community[[from]]$AxA.id, to)
		Community[[to]]$AxA.id <- c(Community[[to]]$AxA.id, from)

		Community[[from]]$AxA.w <- c(Community[[from]]$AxA.w, w)
		Community[[to]]$AxA.w <- c(Community[[to]]$AxA.w, w)

		Community[[from]]$AxA.t <- c(Community[[from]]$AxA.t, t)
		Community[[to]]$AxA.t <- c(Community[[to]]$AxA.t, t)
	}
	return(Community)
}

get.connections <- function(source, Community) {
	links <- vector("list")
	author <- Community[[source]]
	m <- length(author$AxA.id)

	if (m == 0) {
		return(links)
	}

	for (i in 1:m) {
		link <- list(id = author$AxA.id[i], w = author$AxA.w[i])
		links <- c(links, list(link))
	}
	return(links)
}

get.degree.distro <- function(Community) {
	links <- vector("list")
	n <- length(Community)

	if (n == 0) {
		return(links)
	}

	for (i in 1:n) {
		author <- Community[[i]]
		d <- length(author$AxA.id)
		link <- list(id = author$ID, w = d)
		links <- c(links, list(link))
	}

	return(links)
}

get.common.knowledge <- function(Author.List) {
	keys <- vector("integer")
	nA <- length(Author.List)
	if (nA == 0) {
		return(keys)
	}

	keys <- Author.List[[1]]$K.Set
	if (nA == 1) {
		return(keys)
	}

	for (i in 2:nA) {
		author <- Author.List[[i]]
		keys <- intersect(keys, author$K.Set)
	}

	return(keys)
}

get.knowledge.degrees <- function(Author.List, keys = NULL) {
	links <- vector("list")

	if (length(Author.List) == 0) {
		return(links)
	}

	if (is.null(keys)) {
		keys <- vector("integer")
		for (i in 1:length(Author.List)) {
			author <- Author.List[[i]]
			keys <- union(keys, author$K.Set)
		}

		if (length(keys) == 0) {
			return(links)
		}

		for (i in 1:length(keys)) {
			link <- list(id = keys[i], w = 0)
			links <- c(links, list(link))
		}

	} else {
		for (i in 1:length(keys)) {
			link <- list(id = keys[i], w = 0)
			links <- c(links, list(link))
		}
	}

	for (i in 1:length(keys)) {
		key <- keys[i]
		f <- 0

		for (j in 1:length(Author.List)) {
			author <- Author.List[[j]]
			f <- f + length(which(author$KxK.key1 == key))
			f <- f + length(which(author$KxK.key2 == key))
		}

		for (j in 1:length(links)) {
			if (links[[j]]$id == key) {
				links[[j]]$w <- links[[j]]$w + f
			}
		}
	}

	return(links)
}


set.random.walk.probs <- function(links, walk, k = 1) {
	tot.walk <- 0
	m <- length(links)

	for (i in 1:m) {
		tot.walk <- tot.walk + links[[i]]$w
	}
	
	if (m == 1){
		if(is.element(links[[1]]$id, walk)){
			tot.walk <- 100 * tot.walk
		}
	}

	for (i in 1:m) {
		links[[i]]$w <- (links[[i]]$w/tot.walk)^k
	}

	return(links)
}

pref.attach.choice <- function(degrees, is.node = T, total.w = NULL) {
	# degrees is, e.g.,a vector of list(id = 1, w = 1)
	# higer the degree of a node higher the chance for that node to be picked.

	n <- length(degrees)
	if (is.null(total.w)) {
		total.w = 0
		for (i in 1:n) {
			total.w <- total.w + degrees[[i]]$w
		}
	}
	if (is.node) {
		total.w <- total.w + n
	}

	location <- runif(1, 0, total.w)
	ans <- NULL
	for (i in 1:n) {
		width <- degrees[[i]]$w
		if (is.node) {
			width <- width + 1
		}
		if (location < width) {
			ans <- degrees[[i]]$id
			break
		}
		location <- location - width
	}
	return(ans)
}

get.a.walk <- function(Net, start, k = 0, path = vector("integer")) {
	if (length(start) == 0) {
		return(c(path, start))
	}
	repeat {
		path <- union(path, start)
		links <- get.connections(start, Net)
		if (length(links) == 0) {
			break
		}
		k <- k + 1
		link.probs <- set.random.walk.probs(links, path, k)
		start = pref.attach.choice(link.probs, is.node = F, total.w = 1)
		if (is.null(start)) {
			break
		}
	}
	return(path)
}

get.a.semantic.walk <- function(Author, start, k = 0, path = vector("integer")) {
	
	if (!is.element(start, Author$K.Set)) {
		return(path)
	}

	repeat {
		path <- union(path, start)
		links <- get.knowledge.ego.net(Author, start)
		if (length(links) == 0) {
			break
		}
		k <- k + 1
		link.probs <- set.random.walk.probs(links, path, k)
		start = pref.attach.choice(link.probs, is.node = F, total.w = 1)
		if (is.null(start)) {
			break
		}
	}
	return(path)
}


is.a.link <- function(a, b, k1.list, k2.list) {
	# given lists of links indices, checks if a-b link exists. 
	indices.1 <- intersect(which(k1.list == a), which(k2.list == b))
	indices.2 <- intersect(which(k1.list == b), which(k2.list == a))
	return(union(indices.1, indices.2))
}

add.knowledge <- function(KxK.data, Authors) {
	nK <- length(KxK.data[[1]])

	for (i in 1:nK) {
		a.id <- KxK.data$author[i]
		k1 <- KxK.data$k1[i]
		k2 <- KxK.data$k2[i]
		w <- KxK.data$w[i]
		t <- KxK.data$t[i]

		author <- Authors[[a.id]]
		if (k1 == 0) {
			author$K.Set <- union(author$K.Set, k2)
			Authors[[a.id]] <- author
			next
		}
		if (k2 == 0) {
			author$K.Set <- union(author$K.Set, k1)
			Authors[[a.id]] <- author
			next
		}
		author$K.Set <- union(author$K.Set, c(k1, k2))
		author$KxK.key1 <- c(author$KxK.key1, k1)
		author$KxK.key2 <- c(author$KxK.key2, k2)
		author$KxK.w <- c(author$KxK.w, w)
		author$KxK.t <- c(author$KxK.t, t)
		Authors[[a.id]] <- author
	}

	return(Authors)
}

get.knowledge.ego.net <- function(Author, key) {

	links <- vector("list")

	if (!is.element(key, Author$K.Set)) {
		return(links)
	}

	key1.list <- Author$KxK.key1
	key2.list <- Author$KxK.key2

	if (length(key1.list) == 0) {
		link <- list(id = key, w = 1)
		links <- c(links, list(link))
		return(links)
	}

	indices.key2 <- which(key1.list == key)
	indices.key1 <- which(key2.list == key)

	node.set <- vector("integer")
	m = length(indices.key1)
	for (i in 1:m) {
		if (m == 0) {
			break
		}
		ind <- indices.key1[i]
		key1 <- Author$KxK.key1[ind]
		node.set <- union(node.set, key1)
	}
	m = length(indices.key2)
	for (i in 1:m) {
		if (m == 0) {
			break
		}
		ind <- indices.key2[i]
		key2 <- Author$KxK.key2[ind]
		node.set <- union(node.set, key2)
	}

	m = length(node.set)
	if (m == 0) {
		return(links)
	}
	for (i in 1:m) {
		link <- list(id = node.set[i], w = 0)
		links <- c(links, list(link))
	}

	n = length(node.set)
	m = length(indices.key1)
	for (i in 1:m) {
		if (m == 0) {
			break
		}
		ind <- indices.key1[i]
		key <- Author$KxK.key1[ind]
		for (j in 1:n) {
			if (key == links[[j]]$id) {
				links[[j]]$w <- links[[j]]$w + Author$KxK.w[ind]
				break
			}
		}
	}

	m = length(indices.key2)
	for (i in 1:m) {
		if (m == 0) {
			break
		}
		ind <- indices.key2[i]
		key <- Author$KxK.key2[ind]
		for (j in 1:n) {
			if (key == links[[j]]$id) {
				links[[j]]$w <- links[[j]]$w + Author$KxK.w[ind]
				break
			}
		}
	}

	return(links)
}

update.AxA <- function(Community, coauthors, t, w = 1) {

	if (length(coauthors) == 0) {
		return(Community)
	}

	if (length(coauthors) == 1) {
		from <- coauthors[1]
		Community[[from]]$Pub.Count <- Community[[from]]$Pub.Count + 1
		return(Community)
	}

	nc <- length(coauthors)

	for (i in 1:nc) {
		ind <- coauthors[i]
		Community[[ind]]$Pub.Count <- Community[[ind]]$Pub.Count + 1
	}

	nc <- nc - 1
	for (i in 1:nc) {
		from <- coauthors[1]
		others <- setdiff(coauthors, from)
		for (j in 1:length(others)) {
			to <- others[j]
			Community[[from]]$A.Set <- union(Community[[from]]$A.Set, to)
			Community[[to]]$A.Set <- union(Community[[to]]$A.Set, from)


			Community[[from]]$AxA.id <- c(Community[[from]]$AxA.id, to)
			Community[[to]]$AxA.id <- c(Community[[to]]$AxA.id, from)

			Community[[from]]$AxA.w <- c(Community[[from]]$AxA.w, w)
			Community[[to]]$AxA.w <- c(Community[[to]]$AxA.w, w)

			Community[[from]]$AxA.t <- c(Community[[from]]$AxA.t, t)
			Community[[to]]$AxA.t <- c(Community[[to]]$AxA.t, t)
		}
		coauthors <- others
	}

	return(Community)
}


update.KxK <- function(Community, coauthors, keywords, t, w = 1) {

	if (length(keywords) == 0) {
		return(Community)
	}

	if (length(coauthors) == 0) {
		return(Community)
	}

	nc <- length(coauthors)
	nk <- length(keywords)

	for (i in 1:nc) {
		a.id <- coauthors[i]
		Community[[a.id]]$K.Set <- union(Community[[a.id]]$K.Set, keywords)

		if (nk == 1) {
			next
		}

		m <- nk - 1
		keys <- keywords
		for (k in 1:m) {
			key1 <- keys[1]
			others <- setdiff(keys, key1)
			for (j in 1:length(others)) {
				key2 <- others[j]
				Community[[a.id]]$KxK.key1 <- c(Community[[a.id]]$KxK.key1, key1)
				Community[[a.id]]$KxK.key2 <- c(Community[[a.id]]$KxK.key2, key2)
				Community[[a.id]]$KxK.w <- c(Community[[a.id]]$KxK.w, w)
				Community[[a.id]]$KxK.t <- c(Community[[a.id]]$KxK.t, t)
			}
			keys <- others
		}
	}
	return(Community)
}

compute.Ksimilarity <- function(A1, A2) {
	sim <- 0
	x <- length(intersect(A1$K.Set, A2$K.Set))
	a <- length(A1$K.Set)
	b <- length(A1$KxK.key1)

	if (a == 0) {
		return(sim)
	}

	if ((length(A2$KxK.key1) == 0) | (b == 0)) {
		similarity <- x/a
		return(sim)
	}

	y <- 0

	for (i in 1:length(A1$KxK.key1)) {
		ak1 <- A1$KxK.key1[i]
		ak2 <- A1$KxK.key2[i]
		for (j in 1:length(A2$KxK.key1)) {
			bk1 <- A2$KxK.key1[j]
			bk2 <- A2$KxK.key2[j]
			if ((ak1 == bk1 & ak2 == bk2) | (ak1 == bk2 & ak2 == bk1)) {
				y <- y + 1
			}
		}
	}

	sim <- ALFA * (x/a) + (1 - ALFA) * (y/b)

	return(sim)
}


compute.Kdiff <- function(A1, A2) {
	diff <- 0
	x <- length(setdiff(A2$K.Set, A1$K.Set))
	a <- length(A2$K.Set)
	b <- length(A2$KxK.key1)

	if (a == 0) {
		return(diff)
	}
	if ((length(A1$KxK.key1) == 0) | (b == 0)) {
		diff <- x/a
		return(diff)
	}

	y <- 0

	for (i in 1:length(A1$KxK.key1)) {
		ak1 <- A1$KxK.key1[i]
		ak2 <- A1$KxK.key2[i]
		for (j in 1:length(A2$KxK.key1)) {
			bk1 <- A2$KxK.key1[j]
			bk2 <- A2$KxK.key2[j]
			if ((ak1 == bk1 & ak2 == bk2) | (ak1 == bk2 & ak2 == bk1)) {
				y <- y + 1
			}
		}
	}

	diff <- ALFA * (x/a) + (1 - ALFA) * (1 - y/b)

	return(diff)
}


form.network <- function(Authors, type = "A") {
	from <- vector("integer")
	to <- vector("integer")

	N <- length(Authors)
	if (N == 0) {
		return(NULL)
	}

	if (type == "A") {
		for (i in 1:N) {
			author <- Authors[[i]]
			f <- author$ID
			m <- length(author$AxA.id)
			if (m == 0) {
				next
			}
			for (j in 1:m) {
				t <- author$AxA.id[j]
				from <- c(from, f)
				to <- c(to, t)
			}
		}

		if (length(to) == 0) {
			return(NULL)
		}
		edgelist <- data.frame(from, to)
		return(edgelist)
	}

	for (i in 1:N) {
		author <- Authors[[i]]
		m <- length(author$KxK.key1)
		if (m == 0) {
			next
		}
		for (j in 1:m) {
			f <- author$KxK.key1[j]
			t <- author$KxK.key2[j]
			from <- c(from, f)
			to <- c(to, t)
		}
	}

	if (length(to) == 0) {
		return(NULL)
	}

	edgelist <- data.frame(from, to)
	return(edgelist)
}



plot.as.network <- function(Edges, n.sizes = NULL, n.labels = NULL, f.name = NULL, v.coords = NULL, type = "A") {

	Net <- as.network(Edges)
	mat <- symmetrize(Net, return.as.edgelist = T)
	a.net <- as.network(mat, directed = FALSE)

	if (type == "A") {
		n.colors <- "red"
	} else {
		n.colors <- "blue"
	}

	if (is.null(n.sizes)) {
		degrees <- degree(a.net, gmode = "graph")
		n.sizes <- (degrees + 1)^0.4
	}

	if (is.null(n.labels)) {
		n.labels <- network.vertex.names(a.net)
	}

	if (is.null(f.name)) {
		dev.new()
		if (length(n.sizes) < 2) {
			if (is.null(v.coords)) {
				v.coords <- gplot(a.net, vertex.cex = n.sizes, label.cex = 0.7, vertex.col = n.colors, gmode = "graph")
			} else {
				gplot(a.net, vertex.cex = n.sizes, label.cex = 0.7, vertex.col = n.colors, coord = v.coords, 
					gmode = "graph")
			}
			return(v.coords)
		}
		if (is.null(v.coords)) {
			v.coords <- gplot(a.net, vertex.cex = n.sizes, label = n.labels, label.cex = 0.7, vertex.col = n.colors, 
				gmode = "graph")
		} else {
			gplot(a.net, vertex.cex = n.sizes, label = n.labels, label.cex = 0.7, vertex.col = n.colors, 
				coord = v.coords, gmode = "graph")
		}
		return(v.coords)
	}

	pngfile <- paste0(PLOTS_DIR, f.name, ".png")
	png(pngfile, width = 1000, height = 1000, pointsize = 18)
	if (length(n.sizes) < 2) {
		if (is.null(v.coords)) {
			v.coords <- gplot(a.net, vertex.cex = n.sizes, label.cex = 0.7, vertex.col = n.colors, gmode = "graph")
		} else {
			gplot(a.net, vertex.cex = n.sizes, label.cex = 0.7, vertex.col = n.colors, coord = v.coords, 
				gmode = "graph")
		}
		dev.off()
		return(v.coords)
	}
	if (is.null(v.coords)) {
		v.coords <- gplot(a.net, vertex.cex = n.sizes, label = n.labels, label.cex = 0.7, vertex.col = n.colors, 
			gmode = "graph")
	} else {
		gplot(a.net, vertex.cex = n.sizes, label = n.labels, label.cex = 0.7, vertex.col = n.colors, coord = v.coords, 
			gmode = "graph")
	}
	dev.off()
	return(v.coords)
}