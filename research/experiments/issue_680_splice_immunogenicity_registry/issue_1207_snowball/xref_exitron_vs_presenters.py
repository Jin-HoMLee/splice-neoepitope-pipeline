import os, io, csv, boto3
from botocore.config import Config

# --- exitron MS-presented neoantigens (Wang et al. Mol Cell 2021, Tables S4+S5) ---
S4_classI = ["MTRAAWPLM","RAAWPLMTR","HGRAESAPR","MAWRMTRAA","RMTRAAWPL","MAVPLLRRM",
 "WMMRKTWMK","SLKVECMPK","MTMDMGRFK","QSSWPCTFK","FQCAFLHRL","FYFGGPQYL","RAAAAAGVR",
 "AAAAGVRRR","AAAAAGVRR","KMRAAAAAG","EMAPSWVLK","KEMAPSWVL","ERRRSLNSF","SAFPFPFER",
 "SFCPSAFPF","SSRSWMLVR","MRCWRTTRM","ATPYSFLSL","LMQATPYSF","MQATPYSFL","IQTPPHDTL","TMWKNHYQK"]
S5_classI = ["WGRLLSEY","VTSMKKPW","MEMKMRKL","KPEQDGDKQ","QQKVPRSA"]
S5_classII = ["TMVSSSPTRLSGLWMT","VARLPLCRREPEP","LVEDCLPNGGAQ","MKANPALTMVSSSPTRL",
 "KFSILTKHKAFNRS","ILTKHKAYNTFSIL","KFSILTKHKAYNT","LTKHKKIHTGETPYK","DTSSSSTGHATPL",
 "KFSILTKHKSFSTF","LHPHPQGAEAAAGLLQRAR","RGEGSARRQAGQGPDSRPG","GRGLHPHPQGAEAAAGLL",
 "LTKHKAYNTFSIL","NPITVRNVGKPFGTAQ","NPALTMVSSSPTRLSGL","ANPALTMVSSSPTRLSGL",
 "PALTMVSSSPTRLSGLW","KPFGTAQALLNIRD","KFSILTKHKAYNTF","KIMKANPALTMV"]
classI = S4_classI + S5_classI          # 33 comparable to our class-I presenters
allI = classI + S5_classII              # 54 total

# --- pull presenter peptides from R2 ---
ep=os.environ["R2_ENDPOINT"]; bk=os.environ["R2_BUCKET"]
s3=boto3.client("s3", endpoint_url=ep, aws_access_key_id=os.environ["R2_ACCESS_KEY_ID"],
    aws_secret_access_key=os.environ["R2_SECRET_ACCESS_KEY"],
    config=Config(signature_version="s3v4", region_name="auto"))

presenters=[]   # (peptide, patient, best_allele, presentation_class, genotype_score)
for pt in ["patient_001","patient_002"]:
    key=f"results/{pt}/predictions/mhc_presentation.tsv"
    obj=s3.get_object(Bucket=bk, Key=key)
    txt=obj["Body"].read().decode("utf-8")
    rd=csv.DictReader(io.StringIO(txt), delimiter="\t")
    n=0
    for row in rd:
        pep=(row.get("peptide") or "").strip().upper()
        if pep:
            presenters.append((pep, pt, row.get("best_allele",""),
                               row.get("presentation_class",""), row.get("genotype_presentation_score","")))
            n+=1
    print(f"  {pt}: {n} presenter rows")

pep_set={p[0] for p in presenters}
print(f"TOTAL presenter peptides: {len(presenters)} rows, {len(pep_set)} unique")
print(f"Exitron peptides: {len(classI)} class-I, {len(allI)} total")
print()

# --- 1. exact match: exitron peptide == a presenter peptide ---
exact=[e for e in allI if e in pep_set]
print("=== 1. EXACT matches (exitron == presenter) ===")
print("  ", exact if exact else "NONE")
print()

# --- 2. exitron class-I peptide is a SUBSTRING of some presenter (unlikely; presenters are short) ---
# --- 3. a presenter peptide is a SUBSTRING of an exitron peptide (class-I core inside class-II region) ---
sub_hits=[]
for e in allI:
    for pep,pt,al,cl,gs in presenters:
        if pep!=e and (pep in e or e in pep):
            sub_hits.append((e, pep, pt, al, cl, gs, "presenter_in_exitron" if pep in e else "exitron_in_presenter"))
print("=== 2/3. SUBSTRING overlaps (either direction) ===")
if sub_hits:
    for e,pep,pt,al,cl,gs,d in sub_hits:
        print(f"   exitron={e}  presenter={pep} ({pt}, {al}, class={cl}, gscore={gs})  [{d}]")
else:
    print("   NONE")
print()

# --- 4. shared gene sanity: do any presenter peptides come from the same exitron genes? (needs gene col) ---
print("=== presenter columns available ===")
obj=s3.get_object(Bucket=bk, Key="results/patient_001/predictions/mhc_presentation.tsv")
hdr=obj["Body"].read().decode("utf-8").splitlines()[0]
print("  ", hdr)
