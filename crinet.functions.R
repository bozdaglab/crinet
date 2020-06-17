# Load libraries
require(data.table)
require(igraph)
require(xlsx)
require(bnlearn)
require(parallel)
require(stats)
require(clusterProfiler)
require(org.Hs.eg.db)
require(R.matlab)
require(tidyr)
require(reshape2)

# Function Declarations
generate_ERs = function(RNA, miRNA, interactions, singleErExp.thr, core.no){
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("interactions", "miRNA"), envir=environment())
  res = pbapply::pblapply(cl = cluster, X = 1:length(unique(interactions$target)), FUN = function(i){
    mr = unique(interactions$target)[i]
    reg = interactions$miRNA[interactions$target %in% mr]
    a = miRNA[1,,drop = F]
    a = a[-1,,drop = F]
    for(j in 1:length(reg)){
      r = reg[j]
      s = interactions$normalizedScores[interactions$target %in% mr & interactions$miRNA %in% r]
      a = rbind(a, miRNA[r,,drop = F]*s)
    }
    sumScore = sum(interactions$normalizedScores[interactions$target %in% mr])
    er = colSums(a)/sumScore
    er = as.data.frame(er)
    return(er)
  })
  parallel::stopCluster(cluster)
  ERs = data.frame(name = unique(interactions$target))
  ERs$ER = res
  rownames(ERs) = unique(interactions$target)
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("ERs", "RNA"), envir=environment())
  res = pbapply::pblapply(cl = cluster, X = 1:nrow(ERs), FUN = function(i){
    mr = rownames(ERs)[i]
    a = ERs$ER[i][[1]]
    b = RNA[mr,, drop = T]
    samples = intersect(rownames(a[a$er!=0,,drop = F]), names(b[b != 0]))
    
    if(length(samples) >= 50){
      a = log2(a[samples,]+0.000001)
      b = log2(RNA[mr,samples]+0.000001)
      corr = cor(a, b)
    }else{
      corr = 0
    }
    return(corr)
  })
  parallel::stopCluster(cluster)
  ERs$singleErExp = unlist(res)
  summary(ERs$singleErExp)
  ERs = ERs[ERs$singleErExp < singleErExp.thr,]
  
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("RNA", "ERs", "singleErExp.thr"), envir=environment())
  res = pbapply::pblapply(cl = cluster, X = 1:nrow(ERs), FUN = function(i){
    mr = rownames(ERs)[i]
    a = ERs$ER[i][[1]]
    b = RNA[mr,, drop = T]
    samples = intersect(rownames(a[a$er!=0,,drop = F]), names(b[b != 0]))
    if(length(samples) >= 50){
      fix.samples = samples
      bootstrapCor = sapply(1:1000, function(bootstrap_iter){
        samples = sample(fix.samples, size = length(fix.samples), replace = T)
        a = log2(a[samples,]+0.000001)
        b = log2(RNA[mr,samples]+0.000001)
        corr = cor(a, b)
        return(corr)
      })
      return(sum(bootstrapCor >= singleErExp.thr)/1000)
    }else{
      return(0)
    }
  })
  parallel::stopCluster(cluster)
  ERs$singleErExp.bootP = unlist(res)
  summary(ERs$singleErExp.bootP)
  ERs = ERs[ERs$singleErExp.bootP < 0.01,]
  return(ERs)
}
generate_candidate_pairs = function(ERs){
  candidate.RNAs = rownames(ERs)
  pairs = as.data.table(t(combn(candidate.RNAs,2)))
  colnames(pairs) = c('RNAi', 'RNAj')
  return(pairs)
}
eliminate_wrt_commonRegulator = function(pairs, interactions, core.no){
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("pairs", "interactions"), envir=environment())
  common.miRs = pbapply::pblapply(cl = cluster, X = 1:nrow(pairs), FUN = function(i){
    RNA1 = pairs$RNAi[i]
    RNA2 = pairs$RNAj[i]
    set1 = unique(interactions$miRNA[interactions$target %in% RNA1])
    set2 = unique(interactions$miRNA[interactions$target %in% RNA2])
    overlap = length(intersect(set1, set2))
    return(overlap)
  })
  parallel::stopCluster(cluster)
  return(unlist(common.miRs))
}
eliminate_wrt_hypergeo = function(pairs, interactions, core.no){
  require(stats)
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("pairs", "interactions"), envir=environment())
  invisible(parallel::clusterEvalQ(cl = cluster, expr = library(stats)))
  hypergeometric = pbapply::pblapply(cl = cluster, X = 1:nrow(pairs), FUN = function(i){
    universe = length(unique(interactions$miRNA))
    RNA1 = pairs$RNAi[i]
    RNA2 = pairs$RNAj[i]
    set1 = unique(interactions$miRNA[interactions$target %in% RNA1])
    set2 = unique(interactions$miRNA[interactions$target %in% RNA2])
    overlap = length(intersect(set1, set2))
    score = phyper(overlap, length(set1), universe-length(set1), length(set2), lower.tail = FALSE)
    return(score)
  })
  parallel::stopCluster(cluster)
  return(unlist(hypergeometric))
}
eliminate_wrt_partialCor = function(RNA, cna, pairs, core.no){
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("pairs", "RNA", "cna"), envir=environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn)))
  res = pbapply::pblapply(cl = cluster, X = 1:nrow(pairs), FUN = function(i){
    RNAi = pairs$RNAi[i]
    RNAj = pairs$RNAj[i]
    a = RNA[RNAi,]
    b = RNA[RNAj,]
    samples = intersect(names(a[a != 0]), names(b[b != 0]))
    if(length(samples)>50){
      
      rna.df = log2(RNA[c(RNAi, RNAj),samples]+0.000001)
      cna.df = rna.df[-(1:2),, drop = F]
      if(RNAi %in% rownames(cna)){
        b = cna[RNAi,]
        if(length(names(b[b != 0]))>150){
          samples = intersect(samples, names(b[b != 0]))
          cna.df = rbind(cna.df[, samples], cna[RNAi,samples, drop = FALSE])
          rna.df = rna.df[,samples, drop = F]
        }
      }
      if(RNAj %in% rownames(cna)){
        b = cna[RNAj,]
        if(length(names(b[b != 0]))>150){
          samples = intersect(samples, names(b[b != 0]))
          cna.df = rbind(cna.df[, samples], cna[RNAj,samples, drop = FALSE])
          rna.df = rna.df[,samples, drop = F]
        }
      }
      if(nrow(cna.df) == 0){
        df = as.data.frame(scale(cbind(t(rna.df))))
        pcor = cor(df[,1], df[,2])
      }else{
        df = as.data.frame(scale(cbind(t(rna.df), t(cna.df))))
        colnames(df)[1] = paste(colnames(df)[1], "_ge", sep = "")  
        colnames(df)[2] = paste(colnames(df)[2], "_ge", sep = "")  
        ci.test.result = bnlearn::ci.test(x = colnames(df)[1], y = colnames(df)[2], z = colnames(df)[3:ncol(df)],
                                          data = df, test = "cor")
        pcor = ci.test.result$statistic 
      }
    }else{
      pcor = 0
    }
    return(pcor)
  })
  parallel::stopCluster(cluster)
  return(unlist(res))
}
eliminate_wrt_partialCorP = function(RNA, cna, pairs, core.no){
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("pairs", "RNA", "cna"), envir=environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn)))
  res = pbapply::pblapply(cl = cluster, X = 1:nrow(pairs), FUN = function(i){
    RNAi = pairs$RNAi[i]
    RNAj = pairs$RNAj[i]
    a = RNA[RNAi,]
    b = RNA[RNAj,]
    samples = intersect(names(a[a != 0]), names(b[b != 0]))
    if(length(samples)>50){
      rna.df = log2(RNA[c(RNAi, RNAj),samples]+0.000001)
      cna.df = rna.df[-(1:2),, drop = F]
      if(RNAi %in% rownames(cna)){
        b = cna[RNAi,]
        if(length(names(b[b != 0]))>150){
          samples = intersect(samples, names(b[b != 0]))
          cna.df = rbind(cna.df[, samples], cna[RNAi,samples, drop = FALSE])
          rna.df = rna.df[,samples, drop = F]
        }
      }
      if(RNAj %in% rownames(cna)){
        b = cna[RNAj,]
        if(length(names(b[b != 0]))>150){
          samples = intersect(samples, names(b[b != 0]))
          cna.df = rbind(cna.df[, samples], cna[RNAj,samples, drop = FALSE])
          rna.df = rna.df[,samples, drop = F]
        }
      }
      if(nrow(cna.df) == 0){
        df = as.data.frame(scale(cbind(t(rna.df))))
        pcor = cor(df[,1], df[,2])
        pvalue = cor.test(df[,1], df[,2])$p.value
      }else{
        df = as.data.frame(scale(cbind(t(rna.df), t(cna.df))))
        colnames(df)[1] = paste(colnames(df)[1], "_ge", sep = "")  
        colnames(df)[2] = paste(colnames(df)[2], "_ge", sep = "")  
        ci.test.result = bnlearn::ci.test(x = colnames(df)[1], y = colnames(df)[2], z = colnames(df)[3:ncol(df)],
                                          data = df, test = "cor")
        pvalue = ci.test.result$p.value
      }
    }else{
      pvalue = 1
    }
    return(pvalue)
  })
  parallel::stopCluster(cluster)
  return(unlist(res))
}
eliminate_wrt_partialCor_bootstrapP = function(RNA, cna, pairs, core.no){
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("pairs", "RNA", "cna"), envir=environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn)))
  res = pbapply::pblapply(cl = cluster, X = 1:nrow(pairs), FUN = function(i){
    RNAi = pairs$RNAi[i]
    RNAj = pairs$RNAj[i]
    a = RNA[RNAi,]
    b = RNA[RNAj,]
    samples = intersect(names(a[a != 0]), names(b[b != 0]))
    fix.samples = samples
    
    if(length(samples)>50){
      bootstrapPartialCor = sapply(1:100, function(bootstrap_iter){
        samples = sample(fix.samples, size = length(fix.samples), replace = T)
        rna.df = log2(RNA[c(RNAi, RNAj),samples]+0.000001)
        cna.df = rna.df[-(1:2),, drop = F]
        if(RNAi %in% rownames(cna)){
          b = cna[RNAi,]
          if(length(names(b[b != 0]))>150){
            samples = intersect(samples, names(b[b != 0]))
            cna.df = rbind(cna.df[, samples], cna[RNAi,samples, drop = FALSE])
            rna.df = rna.df[,samples, drop = F]
          }
        }
        if(RNAj %in% rownames(cna)){
          b = cna[RNAj,]
          if(length(names(b[b != 0]))>150){
            samples = intersect(samples, names(b[b != 0]))
            cna.df = rbind(cna.df[, samples], cna[RNAj,samples, drop = FALSE])
            rna.df = rna.df[,samples, drop = F]
          }
        }
        if(nrow(cna.df) == 0){
          df = as.data.frame(scale(cbind(t(rna.df))))
          pcor = cor(df[,1], df[,2])
        }else{
          df = as.data.frame(scale(cbind(t(rna.df), t(cna.df))))
          colnames(df)[1] = paste(colnames(df)[1], "_ge", sep = "")  
          colnames(df)[2] = paste(colnames(df)[2], "_ge", sep = "")  
          ci.test.result = bnlearn::ci.test(x = colnames(df)[1], y = colnames(df)[2], z = colnames(df)[3:ncol(df)],
                                            data = df, test = "cor")
          pcor = ci.test.result$statistic 
        }
        return(pcor)
      })
    }else{
      bootstrapPartialCor = 0
    }
    pval = 0.05
    n = floor(length(bootstrapPartialCor)*pval)+1
    return(unname(sort(bootstrapPartialCor)[n]))
    
  })
  parallel::stopCluster(cluster)
  return(unlist(res))
}
eliminate_wrt_collectiveRegulation = function(RNA, ERs, pairs, core.no){
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("pairs", "RNA", "ERs"), envir=environment())
  res = pbapply::pblapply(cl = cluster, X = 1:nrow(pairs), FUN = function(i){
    RNAi = pairs$RNAi[i]
    RNAj = pairs$RNAj[i]
    sum.of.exp = log2(RNA[RNAi,]+RNA[RNAj,]+0.000001)
    a = ERs[RNAi, "ER"][[1]]
    b = ERs[RNAj, "ER"][[1]]
    samples = intersect(rownames(a[a$er!=0,,drop = F]), rownames(b[b$er!=0,,drop = F]))
    sum.of.ER = log2(a[samples,, drop = F]+b[samples,, drop = F]+0.000001)
    a = sum.of.exp
    b = sum.of.ER
    samples = intersect(names(a[a != 0]), rownames(b[b$er!=0,,drop = F]))
    if(length(samples)>50){
      cor = cor(sum.of.exp[samples], as.vector(sum.of.ER[samples,]))
    }else{
      cor = 1      
    }
    return(cor)
  })
  parallel::stopCluster(cluster)
  return(unlist(res))
}
eliminate_wrt_collectiveRegulation_bootstrapP = function(RNA, ERs, pairs, core.no){
  cluster = parallel::makeCluster(core.no)
  parallel::clusterExport(cl = cluster, varlist = c("pairs", "RNA", "ERs"), envir=environment())
  res = pbapply::pblapply(cl = cluster, X = 1:nrow(pairs), FUN = function(i){
    RNAi = pairs$RNAi[i]
    RNAj = pairs$RNAj[i]
    sum.of.exp = log2(RNA[RNAi,]+RNA[RNAj,]+0.000001)
    a = ERs[RNAi, "ER"][[1]]
    b = ERs[RNAj, "ER"][[1]]
    samples = intersect(rownames(a[a$er!=0,,drop = F]), rownames(b[b$er!=0,,drop = F]))
    sum.of.ER = log2(a[samples,, drop = F]+b[samples,, drop = F]+0.000001)
    a = sum.of.exp
    b = sum.of.ER
    samples = intersect(names(a[a != 0]), rownames(b[b$er!=0,,drop = F]))
    if(length(samples)>50){
      fix.samples = samples
      bootstrapCor = sapply(1:100, function(bootstrap_iter){
        samples = sample(fix.samples, size = length(fix.samples), replace = T)
        cor = cor(sum.of.exp[samples], as.vector(sum.of.ER[samples,]))
        return(cor)
      })
    }else{
      bootstrapCor = 1
    }
    pval = 0.01
    n = floor(length(bootstrapCor)*pval)+1
    return(unname(sort(bootstrapCor, decreasing = T)[n]))
  })
  parallel::stopCluster(cluster)
  return(unlist(res))
}
generate_deconv_input = function(pairs, scores, needNormalization, deconv_in){
  pairs$score = scores
  if(needNormalization){
    my.normalize <- function(x){(x-min(x))/(max(x)-min(x))}
    pairs$score = my.normalize(pairs$score)
  }
  df = pairs[, c("RNAi", "RNAj", "score")]
  x2 = df[,c(2,1,3)]
  colnames(x2) = colnames(x2)[c(2,1,3)]
  x = unique(rbind(df, x2))
  df = tidyr::spread(x, RNAj, score, fill = 0,)
  rownames(df) = df$RNAi
  df = df[,-1]
  R.matlab::writeMat(deconv_in, df = as.matrix(df))
  return(df)
}
eliminate_wrt_networkDeconvolutionWithPercentage = function(pairs, deconv_out, df, perc = .5){
  df_out = R.matlab::readMat(deconv_out)
  df_out = as.data.frame(df_out[[1]])
  rownames(df_out) = rownames(df)
  colnames(df_out) = colnames(df)
  a = as.matrix(df_out)
  z = which(a <= 0, arr.ind = T)
  df_new = as.matrix(df_out)
  rownames(df_new) = colnames(df_new) = colnames(df)
  df_new[df_new == 0] = NA
  df_new = melt(df_new, na.rm = T)
  df_new = df_new[, c(2,1,3)]
  colnames(df_new) = c("RNAi", "RNAj", "ranking")
  thr = round(nrow(pairs)*perc)
  df_new2 = head(df_new[order(df_new$ranking, decreasing= T),], n = thr*2)
  remained.candidate.pairs = merge(pairs, df_new2)
  return(remained.candidate.pairs)
}
