## what there is:
# a single coalescent likelihood
# likelihood defined using "infectiousness" C(t)
# data can be provided for one "group"
# data are clipped to match existing data input
# exp(1) prior on scaling parameter
# no likelihood contribution made when I=0 or t<seed time

## what would be good:
# likelihood calculated over multiple putative trees for each "group"
# likelihoods for multiple groups (and/or their parent(s))
# independent starting time(s), allowing for tree data to start before other data
# independent ending time(s), allowing clipping at end of exponential growth
# comprehensive testing 

###################################################################################################

library(epidemia)
library(lubridate)

data("EuropeCovid")

args <- EuropeCovid
args$algorithm <- "sampling"
args$prior <- shifted_gamma(shape = 1 / 6,
                            scale = 1,
                            shift = log(1.05) / 6)
args$rt <- epirt(
  formula=R(country, date) ~ 0 #+ public_events
)

epp <- epim(rt=args$rt,
            obs=args$obs,
            data=subset(args$data,country%in%c('Sweden','Switzerland')),
            pops=subset(args$pops,country%in%c('Sweden','Switzerland')),
            si=args$si,
            algorithm=args$algorithm,
            sampling_args=list(chains=1))#do.call("epim", args)
plot_rt(epp, group = "Sweden")

###################################################################################################

## coalescent functions
.region_clades2 <- function(td, type = c('ACCTRAN', 'MPR')){
  library (ape )
  library( treedater )
  library( phangorn )
  tr2 <- td ; class(tr2) = 'phylo'
  tr2$tip.label <- unname( tr2$tip.label ) 
  
  
  if (! ( 'sts' %in% names(tr2))){
    sts = as.numeric( sapply( strsplit( tr2$tip.label, '\\|'), function(x) tail(x, 2)[1] ) ) 
    tr2$sts = sts 
  }
  if (! ('Ti' %in% names(tr2) ) ){
    nel = node.depth.edgelength( tr2 ) 
    nel <- nel - max(nel) + max( tr2$sts )
    tr2$Ti <-  nel [ (Ntip(tr2)+1) : (Ntip(tr2)+Nnode(tr2) ) ]
  }
  
  region <- sapply( strsplit( tr2$tip.label, '_|\\|'), function(x) tail(x,1))
  names(region) <- tr2$tip.label
  
  region_levels <- c('Il', 'exog')
  region_pd <- as.phyDat(region, type ='USER' , levels = region_levels )
  
  ancestral.pars(tr2, region_pd , type = type[1] ) -> region_ap
  #plotAnc(tr2,  region_ap, cex = 0) 
  
  nodeStates = ns <- t( sapply(  subset(region_ap, select = 1, site.pattern = TRUE), as.numeric ) )
  iedge =  postorder( tr2 )
  
  ismrca_region = rep(FALSE, Nnode(tr2) + Ntip( tr2 ))
  rootnode = tr2$edge[tail(iedge,1),1]
  
  #region_ancestor 
  clade_mrca = rep( FALSE, Nnode(tr2) + Ntip( tr2 ))
  
  if ( nodeStates[rootnode,1]==1 )
    clade_mrca[rootnode] <- TRUE 
  
  #postorder traversal 
  for ( k in 1:length( iedge )){ 
    i <- iedge[k] 
    uv = tr2$edge[i, ]
    u <- uv[1]
    v = uv[2]
    isexog_u <- nodeStates[u,2]==1
    isexog_v <- nodeStates[v,2]==1
    isregion_u <- ( nodeStates[u,1]==1)
    isregion_v <- ( nodeStates[v,1]==1)
    
    if ( (!isregion_u) & isregion_v & (v>Ntip(tr2)))
      clade_mrca[v] <- TRUE 
    
  }
  clade_mrcas <- which( clade_mrca ) 
  region_clades = lapply( clade_mrcas, function(u){
    tr= extract.clade(tr2, u ) 
    drop.tip( tr, tr$tip.label[ grepl(tr$tip.label, patt = 'exog$') ] )
  })
  # deduplicate overlapping clades 
  region_clades2 <- list() 
  for ( k in 1:length(region_clades) ){
    tr <- region_clades[[k]]
    for (l in (1:length(region_clades))[-k] ){
      tr1 <- region_clades[[l]] 
      if ( Ntip ( tr1 ) < Ntip( region_clades[[k]] )){
        todrop <- intersect( tr$tip.label, tr1$tip.label )
        tr <- drop.tip( tr, todrop )
      }
    }
    region_clades2[[k]] <- tr 
  }
  
  Ti = c( tr2$sts, tr2$Ti )
  list( tres = region_clades2, ancestors = clade_mrcas, tmrcas = Ti[ clade_mrcas ] )
}
.bind_tres <- function(tres, sts){
  phys <- lapply( tres, function(tr) {class(tr) <- 'phylo'; tr } )
  rts <- sapply( phys, function( phy ) {
    mdepth <- max( node.depth.edgelength( phy )  ) 
    mst <- max( sts [ phy$tip.label  ] )
    mst - mdepth  
  })
  minrt <- min( rts ) 
  maxrt <- max( rts ) 
  rels <- rts - minrt + 1e-3 
  phys <- lapply( 1:length(phys), function(k) {
    phy <- phys[[k]]
    phy$root.edge <- rels[k]; phy 
  })
  
  .phy <- phys[[1]]
  if(length( phys )>1)
    for ( k in 2:length( phys )){
      #.phy <- bind.tree( .phy , phys[[k]] )
      .phy <- .phy + phys[[k]] 
    }
  multi2di( .phy )
}
tre2df <- function (apephylo, tre, res, maxHeight = Inf, minLTT = 1, adapt_time_axis = TRUE, 
                    sampleTimes = NULL){
  n <- ape::Ntip(apephylo)
  if (is.null(minLTT)) 
    minLTT <- floor(n/5)
  if (inherits(tre, c("multiPhylo", "list"))) {
    phys <- lapply(tre, function(tr) {
      class(tr) <- "phylo"
      tr
    })
    stopifnot(!is.null(sampleTimes))
    rts <- sapply(phys, function(phy) {
      mdepth <- max(node.depth.edgelength(phy))
      mst <- max(sampleTimes[phy$tip.label])
      mst - mdepth
    })
    mst = max(sampleTimes)
    sts <- sampleTimes - min(rts)
    rhs = mst - rts
    rh = max(rhs)
    maxHeight <- min(rh, maxHeight)
    shs = rh - sts
    inhs_list <- lapply(1:length(phys), function(k) {
      phy <- phys[[k]]
      rh = rhs[k]
      ndel <- ape::node.depth.edgelength(phy)
      sort(rh - ndel[(Ntip(phy) + 1):(Ntip(phy) + phy$Nnode)])
    })
    inhs <- sort(do.call(c, inhs_list))
  }
  else {
    stopifnot(inherits(tre, c("phylo", "treedater")))
    D <- ape::node.depth.edgelength(apephylo)
    rh <- max(D[1:n])
    rhs = rh
    sts <- D[1:n]
    maxHeight <- min(rh, maxHeight)
    shs <- rh - sts
    inhs <- sort(rh - D[(n + 1):(n + apephylo$Nnode)])
  }
  u_shs <- unique(shs)
  u_inhs <- unique(inhs)
  nnode <- sum(inhs <= maxHeight)
  ne_haxis <- seq(maxHeight/res, maxHeight * (1 - 1/res), le = res - 
                    1)
  if (adapt_time_axis) {
    ne_haxis <- approx(seq(0, 1, length.out = nnode), inhs[inhs <= 
                                                             maxHeight], xout = seq(1/res, 1 - 1/res, length.out = res - 
                                                                                      1))$y
  }
  dh_ne <- diff(c(0, ne_haxis, maxHeight))
  tredat <- data.frame(h = c(u_shs, u_inhs, ne_haxis), type = c(rep("sample", 
                                                                    length(u_shs)), rep("node", length(u_inhs)), rep("neswitch", 
                                                                                                                     length(ne_haxis))))
  tredat <- tredat[tredat$h <= maxHeight, ]
  tredat$ne_bin <- sapply(tredat$h, function(x) sum(ne_haxis < 
                                                      x) + 1)
  ltt.h <- function(h) max(1, sum(shs < h) - sum(inhs < h) - 
                             sum(rhs < h))
  tredat$ltt <- sapply(tredat$h, ltt.h)
  tredat$nco <- 0
  tredat$nco[tredat$type == "node"] <- sapply(tredat$h[tredat$type == 
                                                         "node"], function(x) sum(x == inhs))
  tredat <- tredat[order(tredat$h), ]
  tredat$intervalLength <- c(0, diff(tredat$h))
  tredat$ltt_terms <- tredat$ltt * (tredat$ltt - 1)/2
  done <- (minLTT == 1)
  while (!done) {
    done <- TRUE
    bin2maxltt <- sapply(1:res, function(bin) {
      i = which(tredat$ne_bin == bin)
      if (length(i) > 0) 
        return(max(tredat$ltt[i]))
      -Inf
    })
    for (i in 1:nrow(tredat)) {
      if (bin2maxltt[tredat$ne_bin[i]] < minLTT) {
        tredat$ne_bin[i] <- max(1, tredat$ne_bin[i] - 
                                  1)
        done <- FALSE
      }
    }
  }
  tredat$dh <- dh_ne[tredat$ne_bin]
  tredat
}
trdat <- function(td,data_date){
  sampleTimes <- td$sts
  tre = .region_clades2( td )$tres
  
  apephylo <- tre
  if (inherits(tre, c("multiPhylo", "list"))) {
    apephylo <- .bind_tres(tre, sampleTimes)
    class(apephylo) <- "phylo"
  }
  class(apephylo) <- "phylo"
  if (!is.null(sampleTimes)) {
    sampleTimes <- sampleTimes[apephylo$tip]
  }
  
  ## to make number of bins match existing variable N2 in date and length
  t1 <- max(unname(sampleTimes))
  t0 <- max(data_date, min(unname(sampleTimes)))
  tree_delay <- max(floor((min(unname(sampleTimes)) - data_date)*366),0)
  tredat <- tre2df(apephylo = apephylo, tre = tre, res = floor((t1-t0)*366), 
                   maxHeight = t1-t0, adapt_time_axis = F, 
                   sampleTimes = sampleTimes)
  tredat$date <- as.Date(date_decimal(t1-tredat$h))
  #tredat$ne_bin <- tredat$ne_bin + tree_delay
  print(tree_delay)
  list(tredat,tree_delay)
}

## test ###########################################################################################

countries <- c('Spain','Sweden','Switzerland','Italy')
path <- 'trees/'
pl <- list()
for(c0 in 1:length(countries)){
  pl[[c0]] <- list()
  c1 <- countries[c0]
  
  ## first: as in the original, without tree
  epp <- epim(rt=args$rt,
              obs=args$obs,
              data=subset(args$data,country%in%c1),
              pops=subset(args$pops,country%in%c1),
              si=args$si,
              algorithm=args$algorithm,
              sampling_args=list(chains=3))
  print(epp)
  ## save plot 1
  pl[[c0]][[1]] <- plot_rt(epp, group = c1)
  
  ## second: use tree data
  
  ## get tree data
  treefile <- paste0(path,c1,'/trees.Rds')
  td <- readRDS(treefile)[[1]]
  
  ## sort dates
  ## to make number of bins match existing variable N2 in date and length
  tree_dates <- as.Date(date_decimal(unname(td$sts))) # h=0 at max(tree_dates)
  cat(paste0('Tree dates for ',c1,': ',min(tree_dates),'--',max(tree_dates),'\n\n'))
  epp_dates <- unique(subset(args$data,country%in%c1)$date)
  cat(paste0('Data dates for ',c1,': ',min(epp_dates),'--',max(epp_dates),'\n\n'))
  dats <- trdat(td,data_date=decimal_date(min(epp_dates)))
  
  dat <- dats[[1]] # subset(dats[[1]],date<'2020-03-20')
  
  ## put into stan input format
  phylo <- list()
  nms <- c('ne_bin','ltt_terms','nco','intervalLength')
  for(nm in nms) phylo[[paste0(nm,'_phy')]] <- dat[[nm]]
  phylo$N_phy <- length(phylo$ne_bin_phy)
  phylo$tree_delay <- dats[[2]]
  
  args2 <-args
  args2$rt$phylo <- phylo
  
  epp <- epim(rt=args2$rt,
              obs=args$obs,
              data=subset(args$data,country%in%c1),
              pops=subset(args$pops,country%in%c1),
              si=args$si,
              algorithm=args$algorithm,
              sampling_args=list(chains=3))
  print(epp)
  ## save plot 2
  pl[[c0]][[2]] <- plot_rt(epp, group = c1)

  # take early / exponential growth subset
  # dat <- subset(dats[[1]],date<'2020-03-20')
  # 
  # ## put into stan input format
  # phylo <- list()
  # nms <- c('ne_bin','ltt_terms','nco','intervalLength')
  # for(nm in nms) phylo[[paste0(nm,'_phy')]] <- dat[[nm]]
  # phylo$N_phy <- length(phylo$ne_bin_phy)
  # phylo$tree_delay <- dats[[2]]
  # 
  # args2 <-args
  # args2$rt$phylo <- phylo
  # 
  # epp <- epim(rt=args2$rt,
  #             obs=args$obs,
  #             data=subset(args$data,country%in%c1),
  #             pops=subset(args$pops,country%in%c1),
  #             si=args$si,
  #             algorithm=args$algorithm,
  #             sampling_args=list(chains=3))
  # print(epp)
  # ## save plot 3
  # pl[[c0]][[3]] <- plot_rt(epp, group = c1)
  
}

## plot
library(grid)
library(gridExtra)
gr <- arrangeGrob(grobs=do.call(c,pl),
                  widths=c(1,1))
pdf('epidemia_phylo.pdf',height=10)
grid.draw(gr)
dev.off()
