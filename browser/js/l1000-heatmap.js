
function drawL1000Circos(error, genes, signatureData, coexpressionData) {
  // Get drug
  console.log(signatureData[0].drug);
  var width = 800;
  var circosHeatmap = new Circos({
        container: '#l1000-'+signatureData[0].drug,
        width: width,
        height: width
    });

  // Genes
  genes = genes.map(function(d) {
    return {
      id: d,
      label: d,
      color: 'white',
      len: 6
    }
  });

  // Signature
  signature = signatureData.map(function(d, i) {
    pos = i % 6;
    return {
      block_id: d.gene_symbol,
      start: pos,
      end: pos+1,
      value: d.value
    }
  });

  // Concentration
  concentration = signatureData.map(function(d, i) {
    pos = i % 6;
    return {
      block_id: d.gene_symbol,
      start: pos,
      end: pos+1,
      value: d.concentration
    }
  });

  // Timepoint
  timepoint = signatureData.map(function(d, i) {
    pos = i % 6;
    return {
      block_id: d.gene_symbol,
      start: pos,
      end: pos+1,
      value: d.timepoint
    }
  });

  // Coexpression
  coexpressionData = coexpressionData.filter(function(d) { return Math.abs(d.value) > 0.6 });
  console.log(coexpressionData);
  coexpression = coexpressionData.map(function(d, i) {
    return {
      source: {
        id: d.source_gene_symbol,
        start: 2,
        end: 4
      },
      target: {
        id: d.target_gene_symbol,
        start: 2,
        end: 4
      },
      value: -d.value
    }
  });

  console.log(coexpression);

  circosHeatmap
  .layout(
    genes,
    {
      innerRadius: width / 2 - 200,
      outerRadius: width / 2 - 170,
      ticks: {display: false},
      labels: {
        position: 'center',
        display: true,
        size: 12,
        color: '#000',
        radialOffset: 10
      }
    }
  )
  .heatmap('signature', signature, {
    innerRadius: 1.01,
    outerRadius: 1.15,
    color: 'RdBu',
    min: -0.5,
    max: 0.5
  })
  .histogram('concentration', concentration, {
    innerRadius: 1.16,
    outerRadius: 1.25,
    logScale: true,
    color: 'Blues',
    min: 0.04,
    max: 10
  })
  .chords('coexpression', coexpression, {
    radius: 0.99,
    color: 'RdBu',
    min: -1,
    max: 1,
    tooltipContent: function(d) {
      return d.source.id+' > '+d.target.id+' | Spearman R = '+d.value.toFixed(2);
    }
  })
  .render()
}
