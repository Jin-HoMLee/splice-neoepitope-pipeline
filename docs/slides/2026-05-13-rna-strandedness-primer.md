---
marp: true
theme: default
paginate: true
header: 'RNA-Seq Strandedness — Visual Primer'
footer: '2026-05-13 · splice-neoepitope-pipeline · context: Issue #279'
style: |
  section { font-size: 22px; }
  h1 { color: #1e3a8a; }
  h2 { color: #374151; }
  code { background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }
  .callout { background: #fef3c7; border-left: 4px solid #f59e0b; padding: 8px 12px; margin-top: 12px; font-size: 18px; }
  .ok { background: #d1fae5; border-left: 4px solid #10b981; padding: 8px 12px; margin-top: 12px; font-size: 18px; }
---

# RNA-Seq Strandedness
## A Visual Primer

**Context:** pipeline strand-aware alignment ([Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279))

Why this deck exists: strandedness terminology drifts between tools (HISAT2 `F`/`R`/`FR`/`RF`, htseq `yes`/`reverse`, Salmon `ISF`/`ISR`, featureCounts `0`/`1`/`2`). Picking the wrong flag silently produces wrong junction strand calls — no error, no warning, just bad data downstream.

<svg viewBox="0 0 700 60" xmlns="http://www.w3.org/2000/svg" style="background:#f9fafb;">
  <text x="350" y="35" font-family="ui-sans-serif,system-ui" font-size="18" text-anchor="middle" fill="#374151" font-style="italic">"R2 maps to which strand?" — the question this deck answers visually.</text>
</svg>

---

## 1 · The dsDNA setup

<svg viewBox="0 0 700 240" xmlns="http://www.w3.org/2000/svg" font-family="ui-monospace,monospace" font-size="16">
  <defs>
    <marker id="arrow-blue" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#3b82f6"/>
    </marker>
    <marker id="arrow-red" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#ef4444"/>
    </marker>
    <marker id="arrow-green" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#10b981"/>
    </marker>
  </defs>
  <!-- Top strand (sense, +) -->
  <text x="40" y="80" fill="#3b82f6">5'</text>
  <line x1="70" y1="75" x2="620" y2="75" stroke="#3b82f6" stroke-width="2" marker-end="url(#arrow-blue)"/>
  <text x="200" y="68" fill="#3b82f6" font-size="14">A T G C A G T A C G G T A T T C A G</text>
  <text x="640" y="80" fill="#3b82f6">3'</text>
  <text x="40" y="100" fill="#3b82f6" font-size="13" font-family="ui-sans-serif,system-ui">SENSE (+, coding) — same sequence as mRNA</text>
  <!-- Hydrogen bonds -->
  <g stroke="#9ca3af" stroke-dasharray="2,2">
    <line x1="90" y1="120" x2="90" y2="140"/>
    <line x1="130" y1="120" x2="130" y2="140"/>
    <line x1="170" y1="120" x2="170" y2="140"/>
    <line x1="210" y1="120" x2="210" y2="140"/>
    <line x1="250" y1="120" x2="250" y2="140"/>
    <line x1="290" y1="120" x2="290" y2="140"/>
    <line x1="330" y1="120" x2="330" y2="140"/>
    <line x1="370" y1="120" x2="370" y2="140"/>
    <line x1="410" y1="120" x2="410" y2="140"/>
    <line x1="450" y1="120" x2="450" y2="140"/>
    <line x1="490" y1="120" x2="490" y2="140"/>
    <line x1="530" y1="120" x2="530" y2="140"/>
    <line x1="570" y1="120" x2="570" y2="140"/>
    <line x1="610" y1="120" x2="610" y2="140"/>
  </g>
  <!-- Bottom strand (antisense, -) -->
  <text x="200" y="172" fill="#ef4444" font-size="14">T A C G T C A T G C C A T A A G T C</text>
  <text x="640" y="160" fill="#ef4444">5'</text>
  <line x1="620" y1="155" x2="70" y2="155" stroke="#ef4444" stroke-width="2" marker-end="url(#arrow-red)"/>
  <text x="40" y="160" fill="#ef4444">3'</text>
  <text x="40" y="195" fill="#ef4444" font-size="13" font-family="ui-sans-serif,system-ui">ANTISENSE (−, template) — read by RNA polymerase</text>
  <!-- RNA pol direction arrow -->
  <line x1="100" y1="220" x2="600" y2="220" stroke="#10b981" stroke-width="2" stroke-dasharray="6,3" marker-end="url(#arrow-green)"/>
  <text x="350" y="215" fill="#10b981" font-size="13" text-anchor="middle" font-family="ui-sans-serif,system-ui">RNA pol travels 3'→5' along template, synthesises mRNA 5'→3'</text>
</svg>

- DNA strands are **anti-parallel** (top 5'→3', bottom 3'→5')
- RNA polymerase reads the **antisense template**; the mRNA produced has the **same sequence as the sense strand**

---

## 2 · mRNA orientation

<svg viewBox="0 0 700 220" xmlns="http://www.w3.org/2000/svg" font-family="ui-monospace,monospace" font-size="16">
  <!-- antisense template (top, faded) -->
  <text x="40" y="50" fill="#ef4444" opacity="0.5">3'</text>
  <line x1="70" y1="46" x2="620" y2="46" stroke="#ef4444" stroke-width="1.5" opacity="0.5"/>
  <text x="200" y="40" fill="#ef4444" font-size="13" opacity="0.6">T A C G T C A T G C C A T A A G T C</text>
  <text x="640" y="50" fill="#ef4444" opacity="0.5">5'</text>
  <text x="40" y="68" fill="#6b7280" font-size="12" font-family="ui-sans-serif,system-ui">antisense template (DNA)</text>
  <!-- transcription arrow down -->
  <line x1="350" y1="80" x2="350" y2="110" stroke="#10b981" stroke-width="2" marker-end="url(#arrow-green)"/>
  <text x="365" y="100" fill="#10b981" font-size="13" font-family="ui-sans-serif,system-ui">transcription</text>
  <!-- mRNA -->
  <text x="40" y="145" fill="#3b82f6">5'-cap</text>
  <line x1="90" y1="140" x2="600" y2="140" stroke="#3b82f6" stroke-width="2.5" marker-end="url(#arrow-blue)"/>
  <text x="200" y="133" fill="#3b82f6" font-size="14">A U G C A G U A C G G U A U U C A G</text>
  <text x="605" y="145" fill="#3b82f6" font-size="13">- (A)n 3'</text>
  <text x="40" y="170" fill="#3b82f6" font-size="13" font-family="ui-sans-serif,system-ui">mRNA — single-stranded, 5'→3', poly-A at 3' end</text>
  <text x="40" y="200" fill="#6b7280" font-size="13" font-family="ui-sans-serif,system-ui">→ mRNA sequence = SENSE-strand sequence (T → U)</text>
</svg>

- mRNA matches the **sense** strand (5'→3'), with **U** in place of T
- "Forward / sense / coding strand / mRNA-strand" all refer to the same direction
- Poly-A tail at 3' end is the anchor 10x bead chemistry uses to capture mRNA

---

## 3 · Library prep — where strand info lives

<svg viewBox="0 0 700 280" xmlns="http://www.w3.org/2000/svg" font-family="ui-sans-serif,system-ui" font-size="14">
  <!-- mRNA -->
  <line x1="80" y1="30" x2="500" y2="30" stroke="#3b82f6" stroke-width="2.5" marker-end="url(#arrow-blue)"/>
  <text x="40" y="34" fill="#3b82f6" font-size="13">mRNA 5'</text>
  <text x="505" y="34" fill="#3b82f6" font-size="13">3' (A)n</text>
  <!-- RT arrow -->
  <line x1="290" y1="50" x2="290" y2="75" stroke="#374151" stroke-width="1.5" marker-end="url(#arrow-blue)"/>
  <text x="305" y="68" fill="#374151" font-size="13">reverse transcription</text>
  <!-- first-strand cDNA (antisense to mRNA) -->
  <text x="40" y="100" fill="#ef4444" font-size="13">first-strand cDNA 3'</text>
  <line x1="500" y1="96" x2="180" y2="96" stroke="#ef4444" stroke-width="2.5" marker-end="url(#arrow-red)"/>
  <text x="505" y="100" fill="#ef4444" font-size="13">5' (T)n</text>
  <text x="160" y="118" fill="#ef4444" font-size="12" font-style="italic">antisense to mRNA</text>
  <!-- Two paths: stranded vs unstranded -->
  <line x1="200" y1="135" x2="120" y2="170" stroke="#10b981" stroke-width="1.5" marker-end="url(#arrow-green)"/>
  <line x1="380" y1="135" x2="460" y2="170" stroke="#9ca3af" stroke-width="1.5" stroke-dasharray="4,3" marker-end="url(#arrow-blue)"/>
  <text x="135" y="160" fill="#10b981" font-size="13">stranded protocol</text>
  <text x="430" y="160" fill="#6b7280" font-size="13">unstranded protocol</text>
  <!-- Stranded box -->
  <rect x="40" y="180" width="280" height="80" fill="#d1fae5" stroke="#10b981" rx="4"/>
  <text x="180" y="205" fill="#065f46" font-size="14" text-anchor="middle" font-weight="bold">Strand info PRESERVED</text>
  <text x="180" y="225" fill="#065f46" font-size="12" text-anchor="middle">dUTP marks 2nd strand → destroyed</text>
  <text x="180" y="243" fill="#065f46" font-size="12" text-anchor="middle">10x: only one orientation amplified</text>
  <!-- Unstranded box -->
  <rect x="380" y="180" width="280" height="80" fill="#fef3c7" stroke="#f59e0b" rx="4"/>
  <text x="520" y="205" fill="#92400e" font-size="14" text-anchor="middle" font-weight="bold">Strand info LOST</text>
  <text x="520" y="225" fill="#92400e" font-size="12" text-anchor="middle">Both strands sequenced equally</text>
  <text x="520" y="243" fill="#92400e" font-size="12" text-anchor="middle">Can't tell sense from antisense</text>
</svg>

- **First-strand cDNA is antisense to mRNA** — this is the universal step
- Whether strand info survives depends on **what the protocol does next** (dUTP destruction, ligation chemistry, 10x bead architecture)

---

## 4 · Paired-end: FR vs RF

<svg viewBox="0 0 700 260" xmlns="http://www.w3.org/2000/svg" font-family="ui-sans-serif,system-ui" font-size="14">
  <!-- Gene / sense strand -->
  <line x1="60" y1="40" x2="640" y2="40" stroke="#9ca3af" stroke-width="2" marker-end="url(#arrow-blue)"/>
  <text x="40" y="44" fill="#6b7280" font-size="12">5'</text>
  <text x="645" y="44" fill="#6b7280" font-size="12">3'</text>
  <text x="350" y="28" fill="#374151" font-size="13" text-anchor="middle" font-style="italic">gene on + strand (sense)</text>

  <!-- FR (forward-stranded) -->
  <text x="40" y="80" fill="#1e3a8a" font-weight="bold">FR (forward-stranded)</text>
  <text x="40" y="98" fill="#6b7280" font-size="12">e.g. Ligation; R1 matches gene direction</text>
  <line x1="160" y1="118" x2="240" y2="118" stroke="#3b82f6" stroke-width="3" marker-end="url(#arrow-blue)"/>
  <text x="200" y="135" fill="#3b82f6" font-size="13" text-anchor="middle">R1 →</text>
  <line x1="440" y1="118" x2="360" y2="118" stroke="#ef4444" stroke-width="3" marker-end="url(#arrow-red)"/>
  <text x="400" y="135" fill="#ef4444" font-size="13" text-anchor="middle">← R2</text>

  <!-- separator -->
  <line x1="60" y1="155" x2="640" y2="155" stroke="#e5e7eb" stroke-width="1"/>

  <!-- RF (reverse-stranded) -->
  <text x="40" y="180" fill="#1e3a8a" font-weight="bold">RF (reverse-stranded)</text>
  <text x="40" y="198" fill="#6b7280" font-size="12">e.g. dUTP / TruSeq Stranded mRNA; R1 opposes gene, R2 matches</text>
  <line x1="240" y1="218" x2="160" y2="218" stroke="#ef4444" stroke-width="3" marker-end="url(#arrow-red)"/>
  <text x="200" y="235" fill="#ef4444" font-size="13" text-anchor="middle">← R1</text>
  <line x1="360" y1="218" x2="440" y2="218" stroke="#3b82f6" stroke-width="3" marker-end="url(#arrow-blue)"/>
  <text x="400" y="235" fill="#3b82f6" font-size="13" text-anchor="middle">R2 →</text>
</svg>

- HISAT2 PE flag: `--rna-strandness FR` or `RF` (two-letter)
- **FR**: R1 → sense; **RF**: R1 → antisense (R2 → sense)
- Salmon: `ISF` ≈ FR, `ISR` ≈ RF · htseq: `yes` ≈ FR, `reverse` ≈ RF

---

## 5 · Single-end: F vs R

<svg viewBox="0 0 700 220" xmlns="http://www.w3.org/2000/svg" font-family="ui-sans-serif,system-ui" font-size="14">
  <!-- Gene / sense strand -->
  <line x1="60" y1="40" x2="640" y2="40" stroke="#9ca3af" stroke-width="2" marker-end="url(#arrow-blue)"/>
  <text x="40" y="44" fill="#6b7280" font-size="12">5'</text>
  <text x="645" y="44" fill="#6b7280" font-size="12">3'</text>
  <text x="350" y="28" fill="#374151" font-size="13" text-anchor="middle" font-style="italic">gene on + strand (sense)</text>

  <!-- F (forward) -->
  <text x="40" y="80" fill="#1e3a8a" font-weight="bold">F — read matches gene direction (sense)</text>
  <line x1="200" y1="105" x2="380" y2="105" stroke="#3b82f6" stroke-width="3" marker-end="url(#arrow-blue)"/>
  <text x="290" y="125" fill="#3b82f6" font-size="13" text-anchor="middle">read →</text>

  <!-- separator -->
  <line x1="60" y1="145" x2="640" y2="145" stroke="#e5e7eb" stroke-width="1"/>

  <!-- R (reverse) -->
  <text x="40" y="170" fill="#1e3a8a" font-weight="bold">R — read is the reverse complement of the transcript (antisense)</text>
  <line x1="380" y1="195" x2="200" y2="195" stroke="#ef4444" stroke-width="3" marker-end="url(#arrow-red)"/>
  <text x="290" y="213" fill="#ef4444" font-size="13" text-anchor="middle">← read</text>
</svg>

<div class="callout">⚠️ <b>SE syntax pitfall:</b> HISAT2 SE takes a <b>single letter</b> (<code>F</code> or <code>R</code>) — not <code>FR</code>/<code>RF</code>. Passing <code>RF</code> on SE input is a frequent copy-paste error from PE docs.</div>

---

## 6 · 10x Chromium 3' Gene Expression

<svg viewBox="0 0 700 280" xmlns="http://www.w3.org/2000/svg" font-family="ui-sans-serif,system-ui" font-size="13">
  <!-- Bead -->
  <circle cx="80" cy="80" r="40" fill="#e0e7ff" stroke="#6366f1" stroke-width="2"/>
  <text x="80" y="85" text-anchor="middle" fill="#4338ca" font-size="13" font-weight="bold">bead</text>
  <!-- Oligo on bead -->
  <line x1="120" y1="80" x2="320" y2="80" stroke="#6366f1" stroke-width="2"/>
  <rect x="125" y="70" width="50" height="20" fill="#a5b4fc"/><text x="150" y="84" text-anchor="middle" font-size="11" fill="#1e1b4b">cell BC</text>
  <rect x="180" y="70" width="40" height="20" fill="#c7d2fe"/><text x="200" y="84" text-anchor="middle" font-size="11" fill="#1e1b4b">UMI</text>
  <rect x="225" y="70" width="90" height="20" fill="#ddd6fe"/><text x="270" y="84" text-anchor="middle" font-size="11" fill="#1e1b4b">poly-(dT) primer</text>
  <!-- mRNA hybridised -->
  <line x1="320" y1="110" x2="620" y2="110" stroke="#3b82f6" stroke-width="2.5" marker-end="url(#arrow-blue)"/>
  <text x="290" y="114" fill="#3b82f6" font-size="12">5'</text>
  <text x="625" y="114" fill="#3b82f6" font-size="12">3'</text>
  <text x="470" y="100" fill="#3b82f6" font-size="13" text-anchor="middle" font-style="italic">mRNA (poly-A captured)</text>
  <!-- Reverse arrow for RT -->
  <line x1="450" y1="135" x2="450" y2="158" stroke="#374151" stroke-width="1.5" marker-end="url(#arrow-blue)"/>
  <text x="465" y="152" fill="#374151" font-size="12">RT → first-strand cDNA</text>
  <!-- After RT: insert structure -->
  <text x="40" y="195" fill="#374151" font-size="13" font-weight="bold">Final insert (after second-strand + adapters):</text>
  <rect x="40" y="210" width="60" height="22" fill="#fef3c7"/><text x="70" y="225" text-anchor="middle" font-size="11">Read1 primer</text>
  <rect x="100" y="210" width="50" height="22" fill="#a5b4fc"/><text x="125" y="225" text-anchor="middle" font-size="11">cell BC</text>
  <rect x="150" y="210" width="40" height="22" fill="#c7d2fe"/><text x="170" y="225" text-anchor="middle" font-size="11">UMI</text>
  <rect x="190" y="210" width="30" height="22" fill="#ddd6fe"/><text x="205" y="225" text-anchor="middle" font-size="10">(dT)</text>
  <rect x="220" y="210" width="280" height="22" fill="#dbeafe"/><text x="360" y="225" text-anchor="middle" font-size="12" fill="#1e3a8a" font-weight="bold">cDNA — SENSE strand</text>
  <rect x="500" y="210" width="60" height="22" fill="#fef3c7"/><text x="530" y="225" text-anchor="middle" font-size="11">Read2 primer</text>
  <!-- Read arrows -->
  <line x1="105" y1="252" x2="190" y2="252" stroke="#374151" stroke-width="2" marker-end="url(#arrow-blue)"/>
  <text x="148" y="268" fill="#374151" font-size="12" text-anchor="middle">R1 = 28 bp (BC+UMI)</text>
  <line x1="500" y1="252" x2="220" y2="252" stroke="#10b981" stroke-width="2" marker-end="url(#arrow-green)"/>
  <text x="360" y="268" fill="#10b981" font-size="12" text-anchor="middle" font-weight="bold">R2 = cDNA (sense → forward-stranded)</text>
</svg>

- **R1** = barcode + UMI only → not aligned to genome
- **R2** = reads back into the cDNA, **matches the sense (coding) strand**
- → For HISAT2 SE: `--rna-strandness F` (forward-stranded) — see [10x Tech Note CG000376](https://assets.ctfassets.net/an68im79xiti/awNZTarmwqmwxKcvDv9wv/b05500661e36290b8ce59689ff889ea8/CG000376_TechNote_Antisense_Intronic_Reads_SingleCellGeneExpression_RevA.pdf), [scg_lib_structs 10xChromium3v3](https://scg-lib-structs.readthedocs.io/en/latest/ge/10xChromium3v3.html)

---

## 7 · HISAT2 stranded alignment → XS tag

<svg viewBox="0 0 700 240" xmlns="http://www.w3.org/2000/svg" font-family="ui-monospace,monospace" font-size="13">
  <!-- Unstranded case -->
  <rect x="20" y="20" width="320" height="200" fill="#fef3c7" stroke="#f59e0b" rx="4"/>
  <text x="180" y="42" text-anchor="middle" font-size="14" font-weight="bold" fill="#92400e" font-family="ui-sans-serif,system-ui">unstranded</text>
  <text x="35" y="68" font-size="11" fill="#374151">hisat2 -x idx -U r2.fq -S out.sam</text>
  <line x1="35" y1="85" x2="325" y2="85" stroke="#d1d5db"/>
  <text x="35" y="108" font-size="11" fill="#374151">read1  flag=0  ...  NH:i:1</text>
  <text x="35" y="128" font-size="11" fill="#374151">read2  flag=0  ...  NH:i:1</text>
  <text x="35" y="152" font-size="11" fill="#ef4444">  ↑ no XS tag on most reads</text>
  <text x="35" y="172" font-size="11" fill="#ef4444">  ↑ XS only guessed for canonical</text>
  <text x="35" y="192" font-size="11" fill="#ef4444">    splice motifs (GT-AG)</text>
  <text x="180" y="212" text-anchor="middle" font-size="12" fill="#92400e" font-family="ui-sans-serif,system-ui" font-style="italic">strand info: heuristic</text>

  <!-- Stranded case -->
  <rect x="360" y="20" width="320" height="200" fill="#d1fae5" stroke="#10b981" rx="4"/>
  <text x="520" y="42" text-anchor="middle" font-size="14" font-weight="bold" fill="#065f46" font-family="ui-sans-serif,system-ui">stranded (--rna-strandness F)</text>
  <text x="375" y="68" font-size="11" fill="#374151">hisat2 -x idx --rna-strandness F \</text>
  <text x="385" y="82" font-size="11" fill="#374151">  -U r2.fq -S out.sam</text>
  <line x1="375" y1="95" x2="665" y2="95" stroke="#d1d5db"/>
  <text x="375" y="118" font-size="11" fill="#374151">read1  flag=0  ...  XS:A:+  NH:i:1</text>
  <text x="375" y="138" font-size="11" fill="#374151">read2  flag=0  ...  XS:A:-  NH:i:1</text>
  <text x="375" y="162" font-size="11" fill="#10b981">  ↑ XS set from read strand</text>
  <text x="375" y="182" font-size="11" fill="#10b981">  ↑ correct for all spliced</text>
  <text x="375" y="202" font-size="11" fill="#10b981">    reads, regardless of motif</text>
</svg>

- `XS:A:+/-` is the **splice strand** tag — required by every downstream junction tool
- Without `--rna-strandness`, HISAT2 sets XS only when it can guess from splice motifs (GT-AG); reads with non-canonical or ambiguous splices get **no XS**

---

## 8 · regtools `-s XS` + Issue #279 fix

<svg viewBox="0 0 700 220" xmlns="http://www.w3.org/2000/svg" font-family="ui-monospace,monospace" font-size="13">
  <text x="20" y="35" font-size="13" fill="#374151" font-family="ui-sans-serif,system-ui"><tspan font-weight="bold">regtools junctions extract</tspan> — what <code>-s XS</code> does:</text>
  <rect x="20" y="50" width="660" height="40" fill="#f3f4f6" stroke="#d1d5db" rx="3"/>
  <text x="30" y="76" font-size="13" fill="#1f2937">regtools junctions extract -s XS -a 8 -m 50 -M 500000 -o out.bed in.bam</text>

  <text x="20" y="115" font-size="13" fill="#374151" font-family="ui-sans-serif,system-ui">→ regtools reads the <code>XS</code> tag on each spliced read to call junction strand.</text>
  <text x="20" y="135" font-size="13" fill="#374151" font-family="ui-sans-serif,system-ui">→ <b>No XS → strand defaults to "?" or heuristic guess.</b></text>

  <rect x="20" y="155" width="660" height="60" fill="#fee2e2" stroke="#ef4444" rx="4"/>
  <text x="350" y="178" text-anchor="middle" font-size="14" fill="#991b1b" font-weight="bold" font-family="ui-sans-serif,system-ui">Pipeline today (unstranded HISAT2) → unreliable XS → ambiguous strand calls</text>
  <text x="350" y="200" text-anchor="middle" font-size="13" fill="#991b1b" font-family="ui-sans-serif,system-ui">Junction strand may flip for non-canonical splices; downstream filters lose specificity.</text>
</svg>

**Fix proposed in [Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279):** pass `--rna-strandness F` for 10x R2 SE samples.

<div class="callout">⚠️ <b>Correction to Issue #279 body:</b> the issue currently proposes <code>--rna-strandness RF</code>, but (a) <code>RF</code> is PE syntax; SE takes a single letter, and (b) 10x R2 is <b>forward-stranded</b> (R2 = sense strand), so the correct flag is <code>F</code>, not <code>R</code>. Follow-up correction comment will be posted on the issue.</div>

---

## Cheat sheet

|                          | HISAT2 SE | HISAT2 PE | htseq    | Salmon | featureCounts |
|--------------------------|-----------|-----------|----------|--------|---------------|
| Unstranded               | (omit)    | (omit)    | `no`     | `IU`   | `0`           |
| Forward-stranded (sense) | `F`       | `FR`      | `yes`    | `ISF`  | `1`           |
| Reverse-stranded (anti)  | `R`       | `RF`      | `reverse`| `ISR`  | `2`           |

**This pipeline's case:** 10x 3' GEX v3, R2-only SE, normal samples → **forward-stranded** (R2 = sense) → HISAT2 `--rna-strandness F` · featureCounts `1` · htseq `yes`.

**Sources:** [HISAT2 manual](http://daehwankimlab.github.io/hisat2/manual/) · [10x Tech Note CG000376](https://assets.ctfassets.net/an68im79xiti/awNZTarmwqmwxKcvDv9wv/b05500661e36290b8ce59689ff889ea8/CG000376_TechNote_Antisense_Intronic_Reads_SingleCellGeneExpression_RevA.pdf) · [scg_lib_structs 10xChromium3v3](https://scg-lib-structs.readthedocs.io/en/latest/ge/10xChromium3v3.html) · [Strandedness review (Bfunc Gen 2020)](https://academic.oup.com/bfg/article/19/5-6/339/5837822)
