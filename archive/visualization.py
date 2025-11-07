import pandas as pd
from pyvis.network import Network
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os

# This code was greatly inspired by Giuseppe Saldi

# =========================
# 1. Load and prepare data
# =========================
df = pd.read_excel('network_neset.xlsx')
df = df[df['combined_confidences'] >= 0.90]
df.columns = df.columns.str.strip()
df['strength'] = df['correlation']

# =========================
# 2. Define genes of interest and colors
# =========================
genes_of_interest = {
    'dan': 'blue', 'danr': 'blue', 'bsh': 'orange', 'svp': 'orange', 'erm': 'orange',
    'ap': 'orange', 'pdm3': 'orange', 'bab2': 'orange', 'dve': 'orange', 'zfh1': 'orange'
}

# =========================
# 3. Filter for GOI-GOI interactions only
# =========================
df_goi = df[df['regulator'].isin(genes_of_interest) & df['target'].isin(genes_of_interest)].copy()
df_goi.to_csv('network_goi_subgraph.tsv', sep='\t', index=False)

# =========================
# 4. FlyBase IDs
# =========================
gene_to_fbgn = {
    'dan': 'FBgn0039286', 'danr': 'FBgn0039283', 'bsh': 'FBgn0000529', 'svp': 'FBgn0003651',
    'erm': 'FBgn0031375', 'ap': 'FBgn0267978', 'pdm3': 'FBgn0261588',
    'bab2': 'FBgn0025525', 'dve': 'FBgn0020307', 'zfh1': 'FBgn0004606',
}

# =========================
# 5. Build network
# =========================
net = Network(height='750px', width='100%', notebook=False, cdn_resources='in_line')
norm = mcolors.Normalize(vmin=df_goi['strength'].min(), vmax=df_goi['strength'].max())
cmap = cm.get_cmap('coolwarm')

for _, row in df_goi.iterrows():
    reg = row['regulator']
    tgt = row['target']
    strength = row['strength']
    beta_sign = row.get('beta.sign.sum', 0)

    reg_color = genes_of_interest[reg]
    tgt_color = genes_of_interest[tgt]
    reg_shape = "triangle"
    tgt_shape = "dot"

    net.add_node(reg, label=reg, color=reg_color, shape=reg_shape, size=20,
                 title=f"FlyBase: {reg}", fbgn=gene_to_fbgn.get(reg, ''))
    net.add_node(tgt, label=tgt, color=tgt_color, shape=tgt_shape, size=20,
                 title=f"FlyBase: {tgt}", fbgn=gene_to_fbgn.get(tgt, ''))

    # Edge coloring by sign
    if beta_sign > 0:
        edge_color = 'green'
        edge_title = f"<b>Activation</b><br>Correlation: {strength:.3f}"
    elif beta_sign < 0:
        edge_color = 'purple'
        edge_title = f"<b>Repression</b><br>Correlation: {strength:.3f}"
    else:
        edge_color = mcolors.to_hex(cmap(norm(strength)))
        edge_title = f"Correlation: {strength:.3f}"

    net.add_edge(reg, tgt, value=abs(strength), title=edge_title, color=edge_color, arrows='to')

net.show_buttons(filter_=['physics'])

# =========================
# 6. HTML: Controls + Legends + FlyBase click handler
# =========================
html_str = net.generate_html(notebook=False)

control_html = """<div style="position: fixed; top: 10px; left: 10px; background: white; z-index: 1000; padding: 10px; border: 1px solid #ccc;">
  <div>
    <label for="fontSizeSlider">Font Size:</label>
    <input type="range" id="fontSizeSlider" min="8" max="32" value="14">
    <span id="fontSizeValue">14</span>
  </div>
  <div style="margin-top: 10px;">
    <label for="filterInput">Filter Nodes (comma-separated):</label>
    <input type="text" id="filterInput" placeholder="e.g., dan,pdm3">
    <button id="filterButton">Filter/Highlight</button>
  </div>
</div>
<div style="position: fixed; bottom: 10px; right: 10px; background: white; z-index: 1000; padding: 12px; border: 1px solid #ccc; font-size: 15px;">
  <strong>Legend:</strong>
  <ul style="list-style-type: none; padding: 0;">
    <li><span style="color:green;">&#9650;</span> Green edge: Activation</li>
    <li><span style="color:purple;">&#9650;</span> Purple edge: Repression</li>
    <li><b>Arrow:</b> From regulator to target</li>
  </ul>
</div>"""

click_handler = """
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
    var slider = document.getElementById("fontSizeSlider");
    var fontSizeValue = document.getElementById("fontSizeValue");
    slider.addEventListener("input", function() {
        fontSizeValue.textContent = slider.value;
        network.setOptions({ nodes: { font: { size: parseInt(slider.value) } } });
    });
    document.getElementById("filterButton").addEventListener("click", function(){
        var filterList = document.getElementById("filterInput").value.split(",").map(x => x.trim());
        var allNodes = nodes.get();
        var updates = [];
        allNodes.forEach(function(node) {
            if (filterList.includes(node.label)) {
                node.color = { background: "yellow", border: "black" };
                updates.push(node);
            }
        });
        nodes.update(updates);
    });
    network.on("click", function(params) {
        if (params.nodes.length > 0) {
            var node = nodes.get(params.nodes[0]);
            if (node.fbgn) {
                window.open(`https://flybase.org/reports/${node.fbgn}.html`, '_blank');
            }
        }
    });
});
</script>
"""

html_str = html_str.replace("</body>", control_html + click_handler + "\n</body>")

# =========================
# 7. Save outputs
# =========================
with open("goi_network_subgraph.html", "w") as f:
    f.write(html_str)

print("Saved 'goi_network_subgraph.html'")
