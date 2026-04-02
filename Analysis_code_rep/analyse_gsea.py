
def analyse_gsea(adata, group_name, uns_key, outfile_key):
    in_adata = adata


    #get DEG df
    deg_df = scanpy.get_genes_group_df(
        adata=in_adata,
        group=group_name,
        key=uns_key)

    #run GSEA pre rank
    pre_rank = gp.prerank(deg_df.loc[:,['names', 'logfoldchanges']], gene_sets=HALLMARK_2024)
    pre_rank.res2d.to_csv(f"{outfile_key}_pre_rank.csv")

    pre_rank_term2 = pre_rank.res2d.Term
    axes = pre_rank.plot(terms=pre_rank_term2)
    return axes

    #run gsea enrich

    deg_sig = deg_df[deg_df['pvals_adj'] < 0.05]
    deg_up = deg_sig[deg_sig['logfoldchanges'] > 0]
    deg_down = deg_sig[deg_sig['logfoldchanges'] < 0]

    print("upregulated genes:",deg_up.shape)
    print("downregulated genes:",deg_down.shape)

    enr_up = gp.enrichr(gene_list = deg_up['names'], gene_sets=HALLMARK_2024)
    enr_up.res2d.to_csv(f"{outfile_key}_enrich_up.csv")

    enr_down = gp.enrichr(gene_list = deg_down['names'], gene_sets=HALLMARK_2024)
    enr_down.res2d.to_csv(f"{outfile_key}_enrich_down.csv")

    #plot gsea enrich

    enr_up.res2d['Direction'] = "Upregulated"
    enr_down.res2d['Direction'] = "Downregulated"

    enr_res = pandas.concat([enr_up.res2d, enr_down.res2d])
    
    fig1, ax1 = plt.subplots(figsize=(7,7))
    ax1 = gp.dotplot(p17_KvsL_enr_res,x='Direction',
        title='HALLMARK 2024',  
        x_order = ["Upregulated","Downregulated"],
        size=3,
        show_ring=True,
        cmap = NbDr.reversed(),
        ax=ax1)

    plt.savefig(f"{outfile_key}_enrich_dotplot.png", dpi=200, bbox_inches = "tight")
    plt.close()

    fig2, ax2 = plt.subplots(figsize=(7,7))
    ax2 = gp.barplot(enr_res, x='Direction',
        title='HALLMARK 2024',  
        color=['b','r'],
        ax=ax2)

    plt.savefig(f"{outfile_key}_enrich_barplot.png", dpi=200, bbox_inches = "tight")
    plt.close()




    