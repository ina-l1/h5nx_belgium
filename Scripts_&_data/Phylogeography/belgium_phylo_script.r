# FLUCART phylogenetic and phylogeographic analysis 

# Libraries
    library(tidyverse)
    library(here)
    library(seqinr)
    library(ape)
    library(ggtree)
    library(treeio)
    library(phytools)
    library(HDInterval)
    library(seraphim)
    library(diagram)
    library(png)
    library(grid)
    library(ggplot2)
    library(sf)
    library(patchwork)
    library(MetBrewer)

# Base directory
base_dir <- here("scripts_data", "phylogeography")

# Data curation

## Read metadata
meta_BE <- read.csv(here(base_dir, "data", "metadata", "meta_BE.csv"), sep = ";", header = TRUE)
meta_gisaid <- read.csv(here(base_dir, "data", "metadata", "meta_gisaid.csv"), sep = ";", header = TRUE)

## Read sequences
seq_A <- seqinr::read.fasta(here(base_dir, "data", "sequences", "EA-2020-A_concatenated.fas"))
seq_C <- seqinr::read.fasta(here(base_dir, "data", "sequences", "EA-2020-C_concatenated.fas"))
seq_AB <- seqinr::read.fasta(here(base_dir, "data", "sequences", "EA-2021-AB_concatenated.fas"))
seq_BB <- seqinr::read.fasta(here(base_dir, "data", "sequences", "EA-2022-BB_concatenated.fas"))
seq_all <- c(seq_A, seq_C, seq_AB, seq_BB)

## Clean sequence names
name_clean_func <- function(seq) {
    names(seq) <- names(seq) %>% 
        trimws() %>% 
        gsub("^_+", "", .) %>% gsub("[()']", "_", .)
    seq
}

seq_A <- name_clean_func(seq_A)
seq_C <- name_clean_func(seq_C)
seq_AB <- name_clean_func(seq_AB)
seq_BB <- name_clean_func(seq_BB)

## meta_BE: Split latitude and longitude
meta_BE <- meta_BE %>% separate(
    location_coordinates,
    into = c("latitude", "longitude"),
    sep = ",\\s*",
    convert = TRUE, fill = "right", remove = FALSE
)

## meta_BE: Change EPI ID for "A/Vulpes_vulpes/Belgium/8660_0016/2022" from EPI_ISL_14493893 to EPI_ISL_14493899 (checked on GISAID)
meta_BE$id_gisaid[meta_BE$id_name == "A/Vulpes_vulpes/Belgium/8660_0016/2022"] <- "EPI_ISL_14493899"

## meta_BE: Code wild birds/poultry/hobby/etc. 
table(meta_BE$species_type)

meta_BE <- meta_BE %>% mutate(
    species_group = case_when(
        grepl("^WBS|^WBI", species_type, ignore.case = TRUE) ~ "wild_birds",
        grepl("poultry|hobby|trader", species_type, ignore.case = TRUE) ~ "poultry",
        grepl("mammals", species_type, ignore.case = TRUE) ~ "mammals",
        grepl("captive", species_type, ignore.case = TRUE) ~ "captive_birds",
        TRUE ~ "other"
    )
) %>% relocate(species_group, .after = species_type)

'Wild birds: as you did
Poultry: poultry + hobby + trader
Mammals: wild and domestic mammals
Captive birds (= zoo animals etc.)'

## EA-2020-C: Main clade (C2)
keep_C2 <- c("A/turkey/Italy/21VIR9073-3/2021|EPI_ISL_14760547",
    "A/Withe-tiled_eagle/Estonia/TA2124126-1_21VIR10433-11/2021|EPI_ISL_7778754",
    "A/pheasant/Italy/21VIR11162/2021|EPI_ISL_14760728",
    "A/Northern_gannet/Spain/3372-12_22VIR10586-6/2022|EPI_ISL_15878548",
    "A/goose/Czech_Republic/22608-3T/2021|EPI_ISL_9603901",
    "A/goose/Czech_Republic/25322-205/2021|EPI_ISL_12150666",
    "A/turkey/Italy/21VIR11007-1/2021|EPI_ISL_14760725",
    "A/domestic_duck/England/007588/2022|EPI_ISL_11406402",
    "A/northern_gannet/France/22P020006/2022|EPI_ISL_19135634",
    "A/turkey/Italy/21VIR8817-1/2021|EPI_ISL_7733648",
    "A/pheasant/Italy/21VIR11165/2021|EPI_ISL_14760730",
    "A/mute_swan/Czech_Republic/785/2022|EPI_ISL_12139999",
    "A/chicken/Czech_Republic/23404-2T/2021|EPI_ISL_7626512",
    "A/mallard/Netherlands/21038796-002/2021|EPI_ISL_6761023",
    "A/northern_gannet/Spain/22VIR12774-1_3004-12/2022|EPI_ISL_16507094",
    "A/Dove/Norway/2022-07-3-1/2022|EPI_ISL_18901671",
    "A/mallard/Italy/21VIR8919-2/2021|EPI_ISL_7733591",
    "A/Chicken/Austria/21146859-002/2021|EPI_ISL_17785725",
    "A/chicken/Czech_Republic/23404-4K/2021|EPI_ISL_7626506",
    "A/mute_swan/Czech_Republic/21312/2021|EPI_ISL_6328409",
    "A/gadwall/Netherlands/21038566-002/2021|EPI_ISL_7952128",
    "A/chicken/Italy/21VIR11508/2021|EPI_ISL_14760880",
    "A/Eurasian_eagle-owl/Germany-BY/AI01616/2022|EPI_ISL_16122150",
    "A/Red_fox/Norway/2022-07-1768-2T/2022|EPI_ISL_18455203",
    "A/peafowl/England/061723/2021|EPI_ISL_13369775",
    "A/turkey/Italy/21VIR10847-1/2021|EPI_ISL_14760702",
    "A/Buteo_buteo/Belgium/00331-0007/2022|EPI_ISL_16747602",
    "A_barnacle_goose_Denmark_25242-1.01_2021|EPI_ISL_18788586",
    "A/northern_gannet/France/22P019752/2022|EPI_ISL_19135617",
    "A/white-stork/Nigeria/VRD22-305-Dal03_23VIR8254-98/2022|EPI_ISL_18462526",
    "A/swan/France/21P012384/2021|EPI_ISL_6590766",
    "A/chicken/Northern_Ireland/AI41_21VIR11351-5/2021|EPI_ISL_8810660",
    "A/White-tailed_eagle/Norway/2022-07-1997-1/2022|EPI_ISL_18901644",
    "A/turkey/Italy/21VIR11511/2021|EPI_ISL_14760891",
    "A/chicken/Italy/21VIR9940-5/2021|EPI_ISL_14760624",
    "A/chicken/Italy/21VIR9904-16/2021|EPI_ISL_14760603",
    "A/laying_hen/Italy/22VIR48-1/2022|EPI_ISL_14761261",
    "A/chicken/Italy/21VIR10239/2021|EPI_ISL_7733595",
    "A/chicken/Italy/21VIR10635-1/2021|EPI_ISL_14760677",
    "A/Northern_gannet/Norway/2022-07-1073-1/2022|EPI_ISL_18901609",
    "A/turkey/Italy/21VIR10524-1/2021|EPI_ISL_14760663",
    "A/great_black-backed_gull/Netherlands/21039609-001-002/2021|EPI_ISL_7952119",
    "A/duck/Italy/22VIR1294-1/2022|EPI_ISL_14761407",
    "A/Barnacle_goose/Scotland/075563/2021|EPI_ISL_13486738",
    "A/greylag_goose/Ireland/032969_22VIR1325-6/2021|EPI_ISL_11259279",
    "A/Gannet/Scotland/084490/2022|EPI_ISL_13969424",
    "A/goose/Netherlands/21038940-002/2021|EPI_ISL_7267248",
    "A/buzzard/Germany-NI/AI02060/2022|EPI_ISL_18006943",
    "A/turkey/Italy/21VIR9608-2/2021|EPI_ISL_14760574",
    "A/chicken/Italy/21VIR11074-1/2021|EPI_ISL_14760711",
    "A/turkey/Italy/21VIR10622-1/2021|EPI_ISL_14760667",
    "A/turkey/Italy/21VIR10398/2021|EPI_ISL_14760653",
    "A/chicken/England/072926/2021|EPI_ISL_13370614",
    "A/duck/Poland/H2014_21RS3290-17/2021|EPI_ISL_18458071",
    "A/mute_swan/Czech_Republic/22477-2/2021|EPI_ISL_9603918",
    "A/turkey/Italy/21VIR9693-3/2021|EPI_ISL_14760592",
    "A/Caspian_gull/Netherlands/1/2022|EPI_ISL_9261749",
    "A/chicken/England/008775/2022|EPI_ISL_13370699",
    "A/goose/Czech_Republic/23458-5T/2021|EPI_ISL_9603769",
    "A/mute_swan/Austria/2990-1_22VIR582-2/2021|EPI_ISL_18458075",
    "A/turkey/Italy/21VIR9652-2/2021|EPI_ISL_7733629",
    "A/chicken/Italy/21VIR10394/2021|EPI_ISL_14760651",
    "A/chicken/Czech_Republic/61-2/2022|EPI_ISL_11327089",
    "A/turkey/Italy/21VIR9606-2/2021|EPI_ISL_14760570",
    "A/Red_Fox/Norway/2022-04-22642-2-4/2022|EPI_ISL_19232149",
    "A/goose/Netherlands/21038799-002/2021|EPI_ISL_6761024",
    "A/otter/Scotland/M100839/2022|EPI_ISL_14864641",
    "A/polecat/Netherlands/22002047-001/2022|EPI_ISL_9616266",
    "A/chicken/Italy/21VIR10987-1/2021|EPI_ISL_14760714",
    "A/broiler/Ireland/033734_22VIR1325-20/2021|EPI_ISL_11259271",
    "A/chicken/Italy/21VIR9133-21/2021|EPI_ISL_14760550",
    "A/chicken/Italy/21VIR10386/2021|EPI_ISL_14760661",
    "A/grey_heron/Croatia/132-b/2021|EPI_ISL_18361460",
    "A/grey_heron/Croatia/132-c/2021|EPI_ISL_7892492",
    "A/Barnacle_goose/Scotland/072140/2021|EPI_ISL_13453562",
    "A/guinea_fowl/Italy/21VIR11383-1/2021|EPI_ISL_14760762",
    "A/turkey/Italy/21VIR9510-1/2021|EPI_ISL_7733635",
    "A/hawk/Poland/MB544/2021|EPI_ISL_18245813",
    "A/chicken/Italy/21VIR10616-1/2021|EPI_ISL_14760670",
    "A/turkey/England/007586/2022|EPI_ISL_13370696",
    "A/chicken/Italy/21VIR10812-11/2021|EPI_ISL_14760696",
    "A/turkey/Italy/21VIR11077-1/2021|EPI_ISL_14760720",
    "A/goose/Czech_Republic/22750/2021|EPI_ISL_9603924",
    "A/chicken/England/057314/2021|EPI_ISL_9029956",
    "A/chicken/Italy/21VIR11572/2021|EPI_ISL_14761020",
    "A/chicken/Italy/21VIR10093/2021|EPI_ISL_14760626",
    "A/duck/Czech_Republic/913/2022|EPI_ISL_12140437",
    "A/chicken/Italy/22VIR1953-3/2022|EPI_ISL_14761485",
    "A/great_skua/Scotland/5/2022|EPI_ISL_13453570",
    "A/chicken/Italy/21VIR10358/2021|EPI_ISL_14760640",
    "A/greylag_goose/Ireland/033062_22VIR1325-7/2021|EPI_ISL_11259280",
    "A/European_Herring_Gull/Sweden/SVA220105SZ0249/FB000036/K-2021|EPI_ISL_10937772",
    "A/Morus_bassanus/Belgium/10811_0001/2022|EPI_ISL_16810069",
    "A/mute_swan/Netherlands/21039008-003/2021|EPI_ISL_6675744",
    "A/chicken/Northern_Ireland/AI48_22VIR586-12/2021|EPI_ISL_18458590",
    "A/pheasant/Czech_Republic/25827-1/2021|EPI_ISL_9603931",
    "A/chicken/Iceland/2022AI02565/2022|EPI_ISL_13246657",
    "A/Northern_gannet/Sweden/SVA220525SZ0402/FB001671/O-2022|EPI_ISL_13338081",
    "A/mute_swan/Northern_Ireland/17269/2021|EPI_ISL_13370415",
    "A/Northern_gannet/Norway/2022-07-1073-4/2022|EPI_ISL_18901611",
    "A/turkey/Italy/21VIR9607-2/2021|EPI_ISL_14760569",
    "A/mute_swan/Netherlands/21038362-001/2021|EPI_ISL_6761006",
    "A/Barnacle_goose/Scotland/109090/2021|EPI_ISL_13486785",
    "A/barnacle_goose/Sweden/SVA220104SZ0295/FB000013/M-2022|EPI_ISL_10916124",
    "A/chicken/Italy/21VIR11238/2021|EPI_ISL_14760744",
    "A/turkey/Poland/H1924-T1/2021|EPI_ISL_6930564",
    "A/Barnacle_goose/Scotland/109105/2021|EPI_ISL_13486796",
    "A/Sparrowhawk/Scotland/022623/2022|EPI_ISL_13486962",
    "A/Canada_goose/Luxembourg/21217776/2021|EPI_ISL_13026104",
    "A/Chicken/Iceland/2023AI04237/2022|EPI_ISL_18001810",
    "A/mute_swan/Austria/21154984-002/2021|EPI_ISL_17785729",
    "A/Canada_goose/Luxembourg/21239614/2021|EPI_ISL_13026106",
    "A/chicken/Netherlands/21038675-001005/2021|EPI_ISL_6328036",
    "A/otter/Scotland/M100366/2022|EPI_ISL_14864642",
    "A/northern_gannet/France/22P019853/2022|EPI_ISL_19135625",
    "A/turkey/Italy/21VIR11396-1/2021|EPI_ISL_14760814",
    "A/Kestrel/Scotland/058512/2021|EPI_ISL_13369761",
    "A/chicken/Netherlands/21040811-001-005/2021|EPI_ISL_7952118",
    "A/greater_white-fronted_goose/Sweden/SVA211217SZ0323/FB004845/2021|EPI_ISL_9648238",
    "A/Great_Black-backed_Gull/4/2022|EPI_ISL_15579544",
    "A/turkey/Italy/21VIR11302/2021|EPI_ISL_14760746",
    "A/common_buzzard/Scotland/294149/2021|EPI_ISL_13370418",
    "A/broiler/Italy/21VIR11800/2021|EPI_ISL_14761165",
    "A/turkey/Italy/21VIR9650-2/2021|EPI_ISL_14760586",
    "A/Northern_gannet/Norway/2022-07-1597-1/2022|EPI_ISL_18901639",
    "A/white-tailed_eagle/Germany-MV/AI05975/2021|EPI_ISL_7753272",
    "A/owl/Italy/21VIR11899-1/2021|EPI_ISL_8882215",
    "A/Red_fox/Norway/2022-07-2244-1/2022|EPI_ISL_18901720",
    "A/turkey/Italy/21VIR9648-2/2021|EPI_ISL_14760584",
    "A/swan/Romania/10678_22VIR2749-7/2022|EPI_ISL_11922821",
    "A/Great_cormorant/Luxembourg/21243185/2021|EPI_ISL_13026237",
    "A/wild_bird/Poland/MW542/2021|EPI_ISL_18245812",
    "A/duck/Poland/H1970_21RS3290-7/2021|EPI_ISL_18458061",
    "A/mute_swan/Ireland/033169_22VIR1325-14/2021|EPI_ISL_11259264",
    "A_common_buzzard_Denmark_00526-1.01_2022|EPI_ISL_18788587",
    "A/turkey/Italy/21VIR9891-21/2021|EPI_ISL_14760614",
    "A/sparrowhawk/England/245233/2021|EPI_ISL_13370618",
    "A/chicken/Italy/21VIR11054-1/2021|EPI_ISL_14760736",
    "A/mute_swan/England/244029/2021|EPI_ISL_13369756",
    "A/pheasant/England/294628/2021|EPI_ISL_13370560",
    "A/laying_hen/Italy/22VIR204-1/2022|EPI_ISL_14761318",
    "A_Northern_gannet_Denmark_06250-1.02_2022|EPI_ISL_18788591",
    "A/turkey/Italy/21VIR11492-4/2021|EPI_ISL_14760956",
    "A/turkey/Italy/21VIR9144-2/2021|EPI_ISL_14760549",
    "A/kestrel/Italy/21VIR10468/2021|EPI_ISL_7733642",
    "A/chicken/Italy/22VIR234-1/2022|EPI_ISL_14761312",
    "A/domestic_duck/England/058612/2021|EPI_ISL_9029959",
    "A/Northern_gannet/Norway/2022-07-1116-1/2022|EPI_ISL_18901613",
    "A/Barnacle_goose/Scotland/072152/2021|EPI_ISL_13453563",
    "A/turkey/Italy/21VIR11163/2021|EPI_ISL_14760740",
    "A/northern_gannet/Sweden/SVA220706SZ0236/FB002106/N-2022|EPI_ISL_14986279",
    "A/Grey_heron/Netherlands/1/2022|EPI_ISL_9261748",
    "A/whooper_swan/Ireland/032960_22VIR1325-4/2021|EPI_ISL_11259277",
    "A/goose/Czech_Republic/23458-2T/2021|EPI_ISL_9603767",
    "A/turkey/England/056764/2021|EPI_ISL_8814195",
    "A/mute_swan/Ireland/037311_22VIR1325-13/2021|EPI_ISL_11260218",
    "A/Sandwich_tern/Sweden/SVA220623SZ0227/FB001883/O-2022|EPI_ISL_14966123",
    "A/buzzard/Italy/21VIR11899-5/2021|EPI_ISL_8882216",
    "A/turkey/Poland/H1910-T3/2021|EPI_ISL_6929958",
    "A/Mallard/Netherlands/1/2022|EPI_ISL_10992813",
    "A/pelican/Russia/932-0121/2021|EPI_ISL_16838580",
    "A/african-jacana/Nigeria/VRD22-305-DA19_23VIR8254-114/2022|EPI_ISL_18462541",
    "A/chicken/Italy/21VIR10724-13/2021|EPI_ISL_14760674",
    "A/turkey/Italy/21VIR9520-1/2021|EPI_ISL_14760578",
    "A/chicken/Italy/21VIR11056-6/2021|EPI_ISL_14760738",
    "A/laying_hen/Italy/21VIR11801-7/2021|EPI_ISL_14761176",
    "A/Red-breasted_goose/Jersey/018413/2022|EPI_ISL_13370896",
    "A/gadwall/Netherlands/22000307-001/2022|EPI_ISL_9616210",
    "A/duck/Poland/H1968_21RS3290-5/2021|EPI_ISL_18458059",
    "A/turkey/Italy/21VIR10617-2/2021|EPI_ISL_14760668",
    "A/barnacle_goose/England/294327/2021|EPI_ISL_13370517",
    "A/chicken/England/007584/2022|EPI_ISL_13370697",
    "A/goose/Netherlands/21039029-001/2021|EPI_ISL_7267252",
    "A/Buteo_buteo/belgium/334_0013/2021|EPI_ISL_14388434",
    "A/chicken/England/000187/2022|EPI_ISL_11406401",
    "A/black-headed_gull/Germany-BY/AI01062/2022|EPI_ISL_16096042",
    "A/chicken/Tyumen/27-31V/2021|EPI_ISL_5463793",
    "A/mute_swan/Austria/2990-2_22VIR582-3/2021|EPI_ISL_18458076",
    "A/broiler/Italy/21VIR11847-1/2021|EPI_ISL_14761243",
    "A/Anser_anser/Spain/810-6_22VIR8632-16/2022|EPI_ISL_15234359",
    "A/turkey/Poland/H1911-N/2021|EPI_ISL_6929970",
    "A/chicken/Italy/21VIR10573-1/2021|EPI_ISL_8882169",
    "A/cormorant/Netherlands/22000816-003/2022|EPI_ISL_9616185",
    "A/turkey/Italy/21VIR9768-8/2021|EPI_ISL_7733638",
    "A/northern_gannet/France/22P020333/2022|EPI_ISL_19135729",
    "A/buzzard/Germany-SH/AI06210/2021|EPI_ISL_7753443",
    "A/heron/Italy/22VIR1562/2022|EPI_ISL_14761250",
    "A/chicken/Italy/21VIR10385/2021|EPI_ISL_14760658",
    "A/common_buzzard/Netherlands/22000306-002/2022|EPI_ISL_9616209",
    "A/chicken/Italy/21VIR11073-1/2021|EPI_ISL_14760718",
    "A/white_stork/Germany-MV/AI05979/2021|EPI_ISL_7753329",
    "A/turkey/Italy/21VIR9649-2/2021|EPI_ISL_7733643",
    "A/mute_swan/Latvia/35510-1_23VIR8835-3/2022|EPI_ISL_18509073",
    "A/mute_swan/Croatia/145/2021|EPI_ISL_8338502",
    "A/turkey/Italy/21VIR11661-7/2021|EPI_ISL_14761152",
    "A/turkey/Italy/21VIR9426-3/2021|EPI_ISL_14760563",
    "A/turkey/Italy/21VIR10857/2021|EPI_ISL_14760685",
    "A/northern_gannet/France/22P019834/2022|EPI_ISL_19135623",
    "A/chicken/Italy/21VIR11571/2021|EPI_ISL_14761009",
    "A/turkey/Italy/21VIR11164/2021|EPI_ISL_14760729",
    "A/turkey/Italy/21VIR9651-2/2021|EPI_ISL_14760585",
    "A/turkey/Italy/21VIR10340/2021|EPI_ISL_7733597",
    "A/domestic_duck/England/040831/2022|EPI_ISL_13370920",
    "A/chicken/Italy/21VIR9858-3/2021|EPI_ISL_14760612",
    "A/chicken/Italy/21VIR10634-1/2021|EPI_ISL_14760676",
    "A/laying_hen/Italy/21VIR11168/2021|EPI_ISL_14760732",
    "A/turkey/Italy/21VIR10395/2021|EPI_ISL_14760650",
    "A/turkey/Italy/21VIR9604-3/2021|EPI_ISL_14760571",
    "A/avian/Italy/21VIR11773-4/2021|EPI_ISL_14761208",
    "A/chicken/Italy/21VIR9508-3/2021|EPI_ISL_14760568",
    "A/Great_black-backed_gull/Norway/2022-07-2560-1/2022|EPI_ISL_18901648",
    "A/chicken/Czech_Republic/23404/2021|EPI_ISL_8515480",
    "A/turkey/Italy/21VIR10471/2021|EPI_ISL_14760660",
    "A/chicken/Italy/21VIR11079-11/2021|EPI_ISL_14760716",
    "A/spur-winged_goose/Nigeria/VRD22-305-Ng06_23VIR8254-101/2022|EPI_ISL_18462529",
    "A/chicken/Italy/21VIR10997-1/2021|EPI_ISL_14760723",
    "A/Barnacle_goose/Scotland/075760/2021|EPI_ISL_13486748",
    "A/Chicken/Sweden/SVA220301SZ0012/FB049029/M-2022|EPI_ISL_10910369",
    "A/pheasant/Czech_Republic/25827-2/2021|EPI_ISL_9603937",
    "A/northern_gannet/Spain/22VIR12774-31_3453-1/2022|EPI_ISL_16507107",
    "A/chicken/Italy/21VIR11078-10/2021|EPI_ISL_14760717",
    "A/chicken/Denmark/24357-11/2021-11-02|EPI_ISL_11798570",
    "A/turkey/Italy/21VIR9766-3/2021|EPI_ISL_14760596",
    "A/chicken/England/009811/2022|EPI_ISL_13370701",
    "A/Barnacle_goose/Scotland/003169/2021|EPI_ISL_13486688",
    "A/mute_swan/Czech_Republic/22684/2021|EPI_ISL_9603922",
    "A/Larus_argentatus/Belgium/10811_0002/2022|EPI_ISL_16810105",
    "A/domestic_duck/England/074574/2021|EPI_ISL_13370616",
    "A/chicken/Austria/2786-1_21VIR11354/2021|EPI_ISL_18458073",
    "A/mute_swan/Poland/MB008-22_22VIR5675-16/2022|EPI_ISL_18245597",
    "A/Red_fox/Norway/2022-04-22642-1-4/2022|EPI_ISL_19232148",
    "A/common_kestrel/Netherlands/22000308-002/2022|EPI_ISL_9304991",
    "A/Northern_gannet/Norway/2022-07-1991-1/2022|EPI_ISL_18901642",
    "A/chicken/Italy/21VIR10470/2021|EPI_ISL_14760654",
    "A/domestic_duck/England/069098/2021|EPI_ISL_13370563",
    "A/chicken/England/071013/2021|EPI_ISL_13370606",
    "A/Black-headed_gull/Netherlands/1/2022|EPI_ISL_9261750",
    "A/domestic_goose/Poland/H1931-T1/2021|EPI_ISL_6931008",
    "A/greylag_goose/Netherlands/21038570-002/2021|EPI_ISL_7267243",
    "A/turkey/Italy/21VIR9475-1/2021|EPI_ISL_14760565",
    "A/buzzard/Germany-NI/AI00935/2022|EPI_ISL_14570650",
    "A/duck/Ireland/036105_22VIR1325-22/2021|EPI_ISL_11259273",
    "A/chicken/Netherlands/21038165-006010/2021|EPI_ISL_6101848",
    "A/mute_swan/Netherlands/21039841-002/2021|EPI_ISL_7952124",
    "A/turkey/Italy/21VIR10929/2021|EPI_ISL_14760697",
    "A/little_owl/Italy/21VIR11382-1/2021|EPI_ISL_8882209",
    "A/common_buzzard/England/300284/2021|EPI_ISL_13369770",
    "A/chicken/Italy/21VIR10808-1/2021|EPI_ISL_14760692",
    "A/northern_gannet/Sweden/SVA220719SZ0404/FB002329/O-2022|EPI_ISL_14986700",
    "A/Common_Buzzard/Sweden/SVA211223SZ0343/FB004885/M-2021|EPI_ISL_10886140",
    "A/turkey/England/009687/2022|EPI_ISL_13370700",
    "A/tawny_owl/England/113433/2022|EPI_ISL_13453569",
    "A/Gull/Norway/2022-07-1996-1/2022|EPI_ISL_18901643",
    "A/turkey/Italy/21VIR10092/2021|EPI_ISL_14760630",
    "A/turkey/Italy/21VIR11586-2/2021|EPI_ISL_8882166",
    "A/northern_gannet/Spain/22VIR12774-33_3497-10/2022|EPI_ISL_16507108",
    "A/Mallard/Netherlands/13/2021|EPI_ISL_9012577",
    "A/Anser_albifrons/Belgium/15465_0010/2021|EPI_ISL_9160206",
    "A/turkey/Italy/21VIR10856/2021|EPI_ISL_14760683",
    "A/chicken/Italy/21VIR9951-25/2021|EPI_ISL_7733603",
    "A/chicken/Italy/21VIR11485/2021|EPI_ISL_14760837",
    "A/chicken/Italy/21VIR11056-1/2021|EPI_ISL_14760737",
    "A/duck/Italy/21VIR9373-2/2021|EPI_ISL_14760562",
    "A/goose/Czech_Republic/23458-1K/2021|EPI_ISL_9603780",
    "A/turkey/Italy/21VIR9512-1/2021|EPI_ISL_7733632",
    "A/Northern_gannet/Norway/2022-07-1143-1/2022|EPI_ISL_18901619",
    "A/quail/Italy/21VIR9515-1/2021|EPI_ISL_14760575",
    "A/chicken/England/053826/2022|EPI_ISL_13370924",
    "A/fox/Ireland/028829_22VIR12757-5/2022|EPI_ISL_20077608",
    "A/Barnacle_goose/Scotland/072168/2021|EPI_ISL_13486729",
    "A/mute_swan/Austria/22000062-001/2022|EPI_ISL_17785731",
    "A/turkey/Italy/21VIR8585-1/2021|EPI_ISL_5524384",
    "A/northern_gannet/France/22P020789/2022|EPI_ISL_19135806",
    "A/greylag_goose/England/108485/2021|EPI_ISL_13369809",
    "A/turkey/Italy/21VIR11504/2021|EPI_ISL_14760999",
    "A/chicken/Poland/H1962_21RS3290-4/2021|EPI_ISL_18458058",
    "A/Northern_Gannet/Netherlands/6/2022|EPI_ISL_15267027",
    "A/chicken/Italy/21VIR9509-1/2021|EPI_ISL_14760566",
    "A/Fox/Netherlands/EMC2/2022|EPI_ISL_12069288",
    "A/goose/Scotland/031099/2022|EPI_ISL_13453564",
    "A/turkey/Italy/21VIR10851/2021|EPI_ISL_8882172",
    "A/barnacle_goose/Netherlands/22000066-001/2022|EPI_ISL_9221509",
    "A/chicken/Italy/21VIR9219-6/2021|EPI_ISL_14760557",
    "A/turkey/Italy/21VIR10970-3/2021|EPI_ISL_14760703",
    "A/mute_swan/Netherlands/21038362-002/2021|EPI_ISL_6761007",
    "A/chicken/Italy/21VIR10850/2021|EPI_ISL_8882168",
    "A/chicken/Italy/21VIR10262/2021|EPI_ISL_14760642",
    "A/domestic_duck/Poland/H1942-N/2021|EPI_ISL_6934175",
    "A/goose/Netherlands/22000194-002/2022|EPI_ISL_9616205",
    "A/chicken/England/002070/2022|EPI_ISL_11406398",
    "A/northern_gannet/Sweden/SVA220623SZ0232/FB001982/O-2022|EPI_ISL_14811896",
    "A/chicken/England/068020/2021|EPI_ISL_13370605",
    "A/turkey/Ireland/033674_22VIR1325-19/2021|EPI_ISL_11259269",
    "A/domestic_duck/England/017166/2022|EPI_ISL_13370707",
    "A/barnacle_goose/Netherlands/22002507-002/2022|EPI_ISL_10454874",
    "A/chicken/Scotland/078873/2022|EPI_ISL_13453560",
    "A/greylag_goose/Netherlands/21039746-002/2021|EPI_ISL_7952122",
    "A/turkey/Italy/21VIR10757-1/2021|EPI_ISL_14760700",
    "A/turkey/Italy/21VIR9476-1/2021|EPI_ISL_14760564",
    "A/crane/Poland/MB528/2021|EPI_ISL_18245811",
    "A/Sanderling/Netherlands/1/2022|EPI_ISL_10993196",
    "A/Mute_swan/Scotland/075821/2021|EPI_ISL_13486881",
    "A/broiler/Italy/21VIR11802-1/2021|EPI_ISL_14761188",
    "A/turkey/Italy/21VIR9143-2/2021|EPI_ISL_7733633",
    "A/chicken/Poland/H1928_21RS3290-2/2021|EPI_ISL_18458055",
    "A/turkey/Italy/21VIR9905-2/2021|EPI_ISL_14760607",
    "A/chicken/Czech_Republic/63/2022|EPI_ISL_11327090",
    "A/turkey/England/055251/2021|EPI_ISL_8814146",
    "A/chicken/Italy/21VIR10637-2/2021|EPI_ISL_14760678",
    "A/turkey/Italy/21VIR9647-2/2021|EPI_ISL_14760588",
    "A/chicken/Norway/P156/2021|EPI_ISL_18455216",
    "A/chicken/Italy/21VIR11076-1/2021|EPI_ISL_14760715",
    "A/chicken/Italy/21VIR10618-19/2021|EPI_ISL_14760671",
    "A/chicken/Italy/21VIR10696-7/2021|EPI_ISL_14760688",
    "A/domestic_duck/England/073865/2022|EPI_ISL_13453559",
    "A/swan/Slovenia/2060_22VIR777-4/2021|EPI_ISL_11007531",
    "A/chicken/Czech_Republic/23589-1T/2021|EPI_ISL_7626526",
    "A/Greylag_goose/England/108622/2021|EPI_ISL_13370516",
    "A/laying_hen/Italy/21VIR11300/2021|EPI_ISL_14760735",
    "A/herring_gull/Ireland/033533_22VIR1325-10/2021|EPI_ISL_11259261",
    "A/turkey/Italy/21VIR11691/2021|EPI_ISL_14761119",
    "A/greylag_goose/Netherlands/21038942-001/2021|EPI_ISL_7267250",
    "A/mallard_duck/England/388009/2022|EPI_ISL_13370702",
    "A/turkey/Italy/21VIR9694-2/2021|EPI_ISL_14760589",
    "A/chicken/Italy/21VIR10620-1/2021|EPI_ISL_14760673",
    "A/Greylag_goose/Luxembourg/21217773/2021|EPI_ISL_13026099",
    "A/goose/Italy/Italy/21VIR10458/2021|EPI_ISL_14760656",
    "A/laying_hen/Italy/21VIR11395-1/2021|EPI_ISL_14760803",
    "A/chicken/Netherlands/21039901-001-005/2021|EPI_ISL_7952116",
    "A/swan/Netherlands/21037791-001/2021|EPI_ISL_6328009",
    "A/Northern_gannet/Norway/2022-07-2085-1/2022|EPI_ISL_18901645",
    "A/laying_hen/Italy/21VIR11394-1/2021|EPI_ISL_14760793",
    "A/turkey/Italy/21VIR10103/2021|EPI_ISL_14760634",
    "A/barnacle_goose/Denmark/24273-1.02/2021-10-26|EPI_ISL_11798568",
    "A/chicken/Czech_Republic/22224-3K/2021|EPI_ISL_7224454",
    "A/chicken/Czech_Republic/913/2022|EPI_ISL_12150658",
    "A/duck/Italy/21VIR10523-1/2021|EPI_ISL_14760662",
    "A/gray_heron/Denmark/24326-1.02/2021-10-28|EPI_ISL_11798577",
    "A/Pigeon/Sweden/SVA220701SZ0238/FB002074/C-2022|EPI_ISL_14811931",
    "A/turkey/Italy/21VIR10852/2021|EPI_ISL_14760684",
    "A/mute_swan/Netherlands/21039192-001/2021|EPI_ISL_7267246",
    "A/chicken/Italy/21VIR10811-28/2021|EPI_ISL_14760695",
    "A/barnacle_goose/England/293809/2021|EPI_ISL_13369780",
    "A/mute_swan/Northern_Ireland/17176/2021|EPI_ISL_13369778",
    "A/turkey/Italy/21VIR9692-16/2021|EPI_ISL_14760590",
    "A/avian/Italy/21VIR9425-2/2021|EPI_ISL_14760559",
    "A/chicken/Italy/21VIR10973-4/2021|EPI_ISL_14760705",
    "A/chicken/Poland/H157_22VIR2515-3/2022|EPI_ISL_11922809",
    "A/gannet/Ireland/026433_22VIR12757-4/2022|EPI_ISL_20077607",
    "A/lesser_white-fronted_goose/Germany-MV/AI05973/2021|EPI_ISL_7753309",
    "A/gannet/Scotland/3/2022|EPI_ISL_13453571",
    "A/magpie/Scotland/012805/2022|EPI_ISL_13370895",
    "A/laying_hen/Italy/21VIR11167/2021|EPI_ISL_14760731",
    "A/turkey/Italy/21VIR10456/2021|EPI_ISL_7733613",
    "A/Barnacle_goose/Scotland/076711/2021|EPI_ISL_13486772",
    "A/turkey/Italy/21VIR10855/2021|EPI_ISL_14760687",
    "A/chicken/Iceland/2022AI02564/2022|EPI_ISL_13246267",
    "A/duck/Ireland/036646_22VIR1325-18/2021|EPI_ISL_11259268",
    "A/turkey/Italy/21VIR11692/2021|EPI_ISL_14761129",
    "A/turkey/Italy/21VIR10854/2021|EPI_ISL_14760686",
    "A/mallard/Scotland/057554/2021|EPI_ISL_13369757",
    "A/turkey/Italy/21VIR9215-1/2021|EPI_ISL_14760556",
    "A/waterfowl/Netherlands/21037914-006/2021|EPI_ISL_6375134",
    "A/chicken/Northern_Ireland/17903/2021|EPI_ISL_13370564",
    "A/swan/Lithuania/1220PG1_22VIR7255-1/2022|EPI_ISL_15022831",
    "A/turkey/Italy/21VIR11887-3/2021|EPI_ISL_8882170",
    "A/turkey/Italy/21VIR9903-2/2021|EPI_ISL_14760601",
    "A/northern_gannet/France/22P019683/2022|EPI_ISL_19135610",
    "A/turkey/Italy/21VIR8926-1/2021|EPI_ISL_14760545",
    "A/turkey/England/068583/2021|EPI_ISL_13370559",
    "A/chicken/Italy/21VIR10382/2021|EPI_ISL_14760648",
    "A/turkey/Italy/21VIR11058/2021|EPI_ISL_14760739",
    "A/chicken/Tyumen/27-39V/2021|EPI_ISL_5463794",
    "A/Greylag_goose/Isle_of_Man/010298/2022|EPI_ISL_13370693",
    "A/turkey/Italy/21VIR9519-1/2021|EPI_ISL_14760577",
    "A/turkey/Italy/21VIR10645-15/2021|EPI_ISL_14760680",
    "A/chicken/England/115126/2022|EPI_ISL_15038820",
    "A/chicken/Italy/21VIR10263/2021|EPI_ISL_14760643",
    "A/Great_black-backed_gull/Norway/2022-07-2853-1/2022|EPI_ISL_18901650",
    "A/goose/Ireland/3869_22VIR2064-2/2022|EPI_ISL_11259284",
    "A/Northern_gannet/Norway/2022-07-1116-2/2022|EPI_ISL_18901614",
    "A/turkey/Italy/21VIR9894-3/2021|EPI_ISL_14760621",
    "A/turkey/Italy/21VIR9618-7/2021|EPI_ISL_7733587",
    "A/chicken/Italy/21VIR9902-2/2021|EPI_ISL_14760598",
    "A/chicken/Italy/21VIR10357/2021|EPI_ISL_14760641",
    "A/mute_swan/Northern_Ireland/17404/2021|EPI_ISL_13370509",
    "A/fulmar/Scotland/118673/2022|EPI_ISL_15586131",
    "A/Northern_Gannet/Netherlands/22016851-002/2022|EPI_ISL_16865548",
    "A/goose/Czech_Republic/22608-2T/2021|EPI_ISL_9603819",
    "A/rooster/Slovenia/2039_22VIR777-1/2021|EPI_ISL_11007528",
    "A/chicken/Italy/21VIR9074-10/2021|EPI_ISL_14760548",
    "A/Barnacle_goose/Netherlands/6/2022|EPI_ISL_9261751",
    "A/red_fox/Netherlands/22002992-003/2022|EPI_ISL_10255581",
    "A/turkey/Italy/21VIR9579-10/2021|EPI_ISL_14760583",
    "A/turkey/Italy/21VIR10697-1/2021|EPI_ISL_14760689",
    "A/seagull/Italy/21VIR11259-10/2021|EPI_ISL_8882213",
    "A/turkey/Italy/21VIR10338/2021|EPI_ISL_14760647",
    "A/turkey/Italy/21VIR9767-3/2021|EPI_ISL_7733624",
    "A/chicken/Italy/21VIR9904-18/2021|EPI_ISL_14760605",
    "A/chicken/Italy/21VIR10271/2021|EPI_ISL_14760644",
    "A/Gannet/Scotland/084624/2022|EPI_ISL_13969425",
    "A/duck/Italy/21VIR10809-1/2021|EPI_ISL_14760693",
    "A/chicken/Poland/H1976_21RS3290-9/2021|EPI_ISL_18458063",
    "A/swan/Slovenia/2041_22VIR777-2/2021|EPI_ISL_11007529",
    "A/northern_gannet/France/22P020480/2022|EPI_ISL_19135754",
    "A/White-tailed_eagle/Norway/2022-07-3463-1/2022|EPI_ISL_18901674",
    "A/mute_swan/Netherlands/21039202-001/2021|EPI_ISL_7267247",
    "A/northern_gannet/France/22P019741/2022|EPI_ISL_19135615",
    "A/duck/Czech_Republic/23589-1T/2021|EPI_ISL_7626514",
    "A/duck/Poland/H1969_21RS3290-6/2021|EPI_ISL_18458060",
    "A/avian/Italy/21VIR11301/2021|EPI_ISL_14760752",
    "A/red_fox/Netherlands/22002191-003/2022|EPI_ISL_9616272",
    "A/red_fox/Netherlands/22002722-004/2022|EPI_ISL_10220428",
    "A/mallard/England/056080/2021|EPI_ISL_13369743",
    "A/white-faced_whistling_duck/Nigeria/VRD22-305-Dal13_23VIR8254-108/2022|EPI_ISL_18462535",
    "A/seagull/Italy/21VIR9432-2/2021|EPI_ISL_7733621",
    "A/chicken/Italy/21VIR9904-17/2021|EPI_ISL_14760604",
    "A/turkey/Italy/21VIR9870-3/2021|EPI_ISL_14760599",
    "A/chicken/Italy/21VIR10265/2021|EPI_ISL_14760646",
    "A/turkey/Italy/21VIR9520/2021|EPI_ISL_11007548",
    "A/chicken/Italy/21VIR11075-2/2021|EPI_ISL_14760719",
    "A/swan/Spain/4087-1_22VIR2142-1/2021|EPI_ISL_11259292",
    "A/barnacle_goose/England/293824/2021|EPI_ISL_13369779",
    "A/turkey/Italy/21VIR9816-1/2021|EPI_ISL_7733598",
    "A/northern_gannet/France/22P019845/2022|EPI_ISL_19135624",
    "A/polecat/Netherlands/22002792-005/2022|EPI_ISL_10220429",
    "A/domestic_duck/Northern_Ireland/17961/2021|EPI_ISL_13370568",
    "A/turkey/Italy/21VIR10699-1/2021|EPI_ISL_14760682",
    "A/mute_swan/Czech_Republic/785/2022|EPI_ISL_9603943",
    "A/turkey/Italy/21VIR9644-2/2021|EPI_ISL_14760580",
    "A/swan/Slovenia/2073_22VIR777-6/2021|EPI_ISL_11007533",
    "A/chicken/Italy/21VIR10971-3/2021|EPI_ISL_14760704",
    "A/chicken/Italy/21VIR10094/2021|EPI_ISL_14760629",
    "A/chicken/Italy/21VIR9580-23/2021|EPI_ISL_14760582",
    "A/chicken/Italy/21VIR9868-8/2021|EPI_ISL_14760597",
    "A/Northern_gannet/Norway/2022-07-1146-1T/2022|EPI_ISL_18455202",
    "A/northern_gannet/France/22P019379/2022|EPI_ISL_19135606",
    "A/turkey/Italy/21VIR9999/2021|EPI_ISL_14760617",
    "A/chicken/Poland/H1983_21RS3290-13/2021|EPI_ISL_18458066",
    "A/Great_black-backed_gull/Norway/2022-07-3018-1/2022|EPI_ISL_18901672",
    "A/mute_swan/Poland/MB550/2021|EPI_ISL_18245814",
    "A/Common_Buzzard/Sweden/SVA220127SZ0301/FB000177/M-2022|EPI_ISL_10942551",
    "A/chicken/Italy/21VIR10619-1/2021|EPI_ISL_14760672",
    "A/turkey/England/004667/2022|EPI_ISL_13370692",
    "A/turkey/Italy/21VIR10459-1/2021|EPI_ISL_14760664",
    "A/northern_gannet/France/22P019685/2022|EPI_ISL_19135611",
    "A/northern_gannet/Sweden/SVA220623SZ0232/FB001984/N-2022|EPI_ISL_14811929",
    "A/black-headed_gull/England/388256/2022|EPI_ISL_11561593",
    "A/chicken/Italy/21VIR10989/2021|EPI_ISL_14760709",
    "A/turkey/Italy/21VIR8925-6/2021|EPI_ISL_14760544",
    "A/black-headed_gull/England/306270/2022|EPI_ISL_11406403",
    "A/turkey/Israel/537/2021|EPI_ISL_12749687",
    "A/Branta_canadensis/Belgium/595_0008/2022|EPI_ISL_14389133",
    "A/Lapwing/England/067574/2021|EPI_ISL_13370892",
    "A/northern_fulmar/Sweden/SVA220628SZ0328/FB002021/O-2022|EPI_ISL_14985601",
    "A/turkey/Italy/21VIR11690/2021|EPI_ISL_14761141",
    "A/Mallard/Netherlands/15/2021|EPI_ISL_9012580",
    "A/pink_footed_goose/Scotland/058488/2021|EPI_ISL_13369760",
    "A/pink_footed_goose/Scotland/003617/2022|EPI_ISL_13782461",
    "A/crow/Ireland/035624_22VIR1325-11/2021|EPI_ISL_11259262",
    "A/Mute_Swan/England/294146/2021|EPI_ISL_13370506",
    "A/chicken/Italy/21VIR10387/2021|EPI_ISL_14760655",
    "A/guinea_fowl/Italy/21VIR10351/2021|EPI_ISL_14760638",
    "A/Canada_goose/Luxembourg/22012198/2022|EPI_ISL_13026245",
    "A/turkey/England/006963/2022|EPI_ISL_13370695",
    "A/turkey/Ireland/035425_22VIR1325-17/2021|EPI_ISL_11259267",
    "A/Nothern_gannet/France/22P019153/2022|EPI_ISL_18445745",
    "A/goose/Czech_Republic/22608-1T/2021|EPI_ISL_9603795",
    "A/turkey/Italy/21VIR11586-1/2021|EPI_ISL_14761086",
    "A/mute_swan/Czech_Republic/25702-2/2021|EPI_ISL_9603927",
    "A/chicken/Italy/21VIR10972-3/2021|EPI_ISL_14760706",
    "A/Great_skua/Iceland/2023AI04249/2022|EPI_ISL_18001837",
    "A/magpie/Ireland/032958_22VIR1325-5/2021|EPI_ISL_11259278",
    "A/Guinea_fowl/Italy/21VIR10231/2021|EPI_ISL_14760633",
    "A/mute_swan/England/387121/2021|EPI_ISL_13370619",
    "A/chicken/Denmark/24357-11/2021|EPI_ISL_18953101",
    "A/Great_black-backed_gull/Norway/2022-07-1594-1/2022|EPI_ISL_18901638",
    "A/hawk/Italy/22VIR428-5/2022|EPI_ISL_14761338",
    "A/barnacle_goose/Sweden/SVA211216SZ0382/FB004836/2021|EPI_ISL_9650220",
    "A/Northern_fulmar/Norway/2022-07-1356-1/2022|EPI_ISL_18901629",
    "A/turkey/Italy/21VIR9646-2/2021|EPI_ISL_14760587",
    "A/turkey/Italy/21VIR10396/2021|EPI_ISL_14760649",
    "A/barnacle_goose/England/292151/2021|EPI_ISL_13369744",
    "A/turkey/Italy/21VIR9210-1/2021|EPI_ISL_7733647",
    "A/swan/Romania/10394_22VIR2749-3/2022|EPI_ISL_11922817",
    "A/peregrine_falcon/Germany-SH/AI01386/2022|EPI_ISL_18006766",
    "A/barnacle_goose/Netherlands/21038248-001/2021|EPI_ISL_6761012",
    "A/duck/Italy/21VIR10447/2021|EPI_ISL_7733584",
    "A/turkey/England/112259/2022|EPI_ISL_15115556",
    "A/chicken/Czech_Republic/23589-4/2021|EPI_ISL_7626533",
    "A/turkey/Italy/21VIR9209-1/2021|EPI_ISL_14760551",
    "A/Great_cormorant/Norway/2023-07-124-1/2023|EPI_ISL_18963373",
    "A/laying_hen/Italy/22VIR1521-2/2022|EPI_ISL_14761439",
    "A/quail/Italy/21VIR9474-1/2021|EPI_ISL_14760561",
    "A/turkey/Denmark/24325-25/2021-10-30|EPI_ISL_11798580",
    "A/Phalacrocorax_carbo/Belgium/1734_0002/2022|EPI_ISL_14390385",
    "A/mute_swan/Czech_Republic/22477-1/2021|EPI_ISL_9603916",
    "A/duck/Italy/21VIR9836-15/2021|EPI_ISL_14760608",
    "A/guinea_fowl-chicken/Northern_Ireland/AI37_21VIR11349-2/2021|EPI_ISL_8810659",
    "A/Nothern_gannet/France/22P019553/2022|EPI_ISL_18445749",
    "A/black-headed_gull/England/386075/2021|EPI_ISL_13369774",
    "A/turkey/Italy/21VIR10810-1/2021|EPI_ISL_14760694",
    "A/mute_swan/Croatia/146/2021|EPI_ISL_8439476",
    "A/turkey/Italy/21VIR9695-4/2021|EPI_ISL_14760591",
    "A/broiler/Italy/21VIR11583-1/2021|EPI_ISL_14761041",
    "A/sanderling/Netherlands/21039528-002/2021|EPI_ISL_7952121",
    "A/Chicken/Sweden/SVA211217SZ0002/FB307562/2021|EPI_ISL_9640026",
    "A/Whooper_swan/Scotland/072171/2021|EPI_ISL_13486975",
    "A/White-tailed_eagle/Iceland/2023AI04257/2022|EPI_ISL_18001861",
    "A/turkey/Italy/21VIR9901-3/2021|EPI_ISL_14760631",
    "A/laying_hen/Italy/21VIR11253-1/2021|EPI_ISL_14760750",
    "A/chicken/Italy/21VIR9371-1/2021|EPI_ISL_14760558",
    "A/chicken/England/061529/2021|EPI_ISL_13369773",
    "A/duck/Poland/H1980_21RS3290-11/2021|EPI_ISL_18458065",
    "A/stork/Spain/234-2_22VIR2142-7/2022|EPI_ISL_11259317",
    "A/greylag_goose/Germany-SH/AI06144/2021|EPI_ISL_7753406",
    "A/domestic_duck/England/100990/2022|EPI_ISL_14497837",
    "A/Red-breasted_goose/Jersey/018414/2022|EPI_ISL_13370897",
    "A/Barnacle_goose/Scotland/003588/2022|EPI_ISL_13486713",
    "A/pheasant/Italy/22VIR1953-1/2022|EPI_ISL_14761462",
    "A/turkey/Italy/21VIR10397/2021|EPI_ISL_14760652",
    "A/broiler/Italy/22VIR54-1/2022|EPI_ISL_14761251",
    "A/chicken/Italy/21VIR9691-2/2021|EPI_ISL_7733630",
    "A/turkey/Italy/21VIR9217-1/2021|EPI_ISL_14760554",
    "A/chicken/England/042391/2022|EPI_ISL_13370922",
    "A/chicken/Italy/21VIR11234/2021|EPI_ISL_14760734",
    "A/chicken/Scotland/070931/2021|EPI_ISL_13370607",
    "A/turkey/Italy/21VIR11237/2021|EPI_ISL_14760743",
    "A/fox/Ireland/3866_22VIR2064-1/2022|EPI_ISL_11259283",
    "A/peregrine_falcon/Netherlands/21038169-002/2021|EPI_ISL_6761011",
    "A/Greylag_goose/Sweden/SVA211118SZ0354/FB004554/2021|EPI_ISL_9650684",
    "A/magpie/Italy/21VIR9487-2/2021|EPI_ISL_7733612",
    "A/Northern_gannet/Norway/2022-07-1370-1/2022|EPI_ISL_18901633",
    "A/chicken/Italy/21VIR10384/2021|EPI_ISL_7733607",
    "A/northern_gannet/Spain/22VIR12774-4_3104-40/2022|EPI_ISL_16524596",
    "A/broiler/Italy/22VIR1892-3/2022|EPI_ISL_14761507",
    "A/goose/Czech_Republic/913/2022|EPI_ISL_12150657",
    "A/common_buzzard/England/245061/2021|EPI_ISL_13370615",
    "A/chicken/Italy/21VIR10000/2021|EPI_ISL_14760616",
    "A/Est/Parnu/henharrier/TA21-26003/2021|EPI_ISL_19028670",
    "A/duck/Poland/H2008_21RS3290-16/2021|EPI_ISL_18458070",
    "A/chicken/Italy/21VIR10928/2021|EPI_ISL_14760698",
    "A/chicken/Italy/21VIR11071-2/2021|EPI_ISL_14760712",
    "A/northern_gannet/France/22P020201/2022|EPI_ISL_19135648",
    "A/Great_black-backed_gull/Iceland/2023AI04254/2022|EPI_ISL_18001853",
    "A/chicken/England/058797/2021|EPI_ISL_13369759",
    "A/turkey/Italy/21VIR9696-2/2021|EPI_ISL_14760593",
    "A/Northern_gannet/Norway/2022-07-1116-3/2022|EPI_ISL_18901615",
    "A/greater_canada_goose/Netherlands/21038971-001/2021|EPI_ISL_7267251",
    "A/chicken/Tyumen/27-40V/2021|EPI_ISL_5463795",
    "A/barnacle_goose/Germany-SH/AI06145/2021|EPI_ISL_7753407",
    "A/turkey/Italy/21VIR10636-17/2021|EPI_ISL_14760679",
    "A/chicken/Italy/21VIR11385-1/2021|EPI_ISL_14760774",
    "A/duck/Italy/21VIR11501/2021|EPI_ISL_14760988",
    "A/chicken/Netherlands/21040808-001-005/2021|EPI_ISL_7952117",
    "A/guinea_fowl/Italy/21VIR11298/2021|EPI_ISL_14760745",
    "A/turkey/Italy/21VIR10930/2021|EPI_ISL_14760699",
    "A/red_fox/Netherlands/22002636-005/2022|EPI_ISL_10220427",
    "A/chicken/Italy/21VIR11774-2/2021|EPI_ISL_14761220",
    "A/mute_swan/Poland/MB551/2021|EPI_ISL_18245815",
    "A/chicken/Italy/21VIR9580-22/2021|EPI_ISL_7733609",
    "A/turkey/Denmark/24325-26/2021|EPI_ISL_19237887",
    "A/turkey/Italy/21VIR9839-1/2021|EPI_ISL_14760609",
    "A/barnacle_goose/Scotland/062003/2021|EPI_ISL_13369765",
    "A/laying_hen/Italy/21VIR11299/2021|EPI_ISL_14760751",
    "A/turkey/England/004737/2022|EPI_ISL_11406399",
    "A/northern_gannet/France/22P020030/2022|EPI_ISL_19135635",
    "A/turkey/Italy/21VIR10621-1/2021|EPI_ISL_14760669",
    "A/swan/Slovenia/13_22VIR777-8/2022|EPI_ISL_11007535",
    "A/turkey/Poland/H1923_21RS3290-1/2021|EPI_ISL_18458056",
    "A/Barnacle_goose/Scotland/075827/2021|EPI_ISL_13486760",
    "A/mute_swan/Netherlands/21038479-002/2021|EPI_ISL_6761020",
    "A/turkey/Italy/21VIR11804-1/2021|EPI_ISL_8882167",
    "A/pelican/Russia/1032-1/2021|EPI_ISL_16838579",
    "A/white-tailed_eagle/Ireland/032034_22VIR1325-23/2021|EPI_ISL_11259274",
    "A/Turkey/Sweden/SVA211227SZ0006/FB312473/2021|EPI_ISL_9640083",
    "A/chicken/Italy/21VIR9765-12/2021|EPI_ISL_7733634",
    "A/northern_gannet/France/22P020664/2022|EPI_ISL_19135788",
    "A/chicken/Austria/2782_21VIR11353/2021|EPI_ISL_18458072",
    "A/White_Stork/Sweden/SVA220310SZ0365/FB000564/M-2022|EPI_ISL_11580158",
    "A/turkey/Italy/21VIR9518-1/2021|EPI_ISL_14760576",
    "A/turkey/Italy/21VIR11440-2/2021|EPI_ISL_14760870",
    "A/buzzard/Germany-BB/AI00489/2022|EPI_ISL_18006931",
    "A/Eurasian_teal/Germany-BY/AI05977/2021|EPI_ISL_7753367",
    "A/chicken/Wales/116603/2022|EPI_ISL_15585890",
    "A/swan/Romania/10324_22VIR2749-2/2022|EPI_ISL_11922816",
    "A/turkey/Italy/21VIR11386-1/2021|EPI_ISL_14760784",
    "A/avian/Italy/21VIR10913-3/2021|EPI_ISL_14760713",
    "A/chicken/Italy/21VIR9830-27/2021|EPI_ISL_14760602",
    "A/turkey/Italy/21VIR9645-2/2021|EPI_ISL_14760581",
    "A/northern_gannet/France/22P020184/2022|EPI_ISL_19135639",
    "A/chicken/Italy/21VIR9998/2021|EPI_ISL_14760620",
    "A/Great_black-backed_gull/Norway/2022-07-1235-2/2022|EPI_ISL_18901623",
    "A/swan/Austria/22004012-001/2022|EPI_ISL_17785726",
    "A/chicken/Scotland/068163/2021|EPI_ISL_13370557",
    "A/turkey/Italy/21VIR11500-3/2021|EPI_ISL_14760977",
    "A/northern_gannet/Germany-SH/AI03496/2022|EPI_ISL_14573001",
    "A/eurasian_sparrowhawk/Netherlands/22000555-002/2022|EPI_ISL_9117232",
    "A/goose/Czech_Republic/22608-2/2021|EPI_ISL_8515479",
    "A/mute_swan/Ireland/032363_22VIR1325-1/2021|EPI_ISL_11259270",
    "A/turkey/Poland/H1941_21RS3290-3/2021|EPI_ISL_18458057",
    "A/Common_buzzard/Scotland/020817/2022|EPI_ISL_13486807",
    "A/northern_gannet/France/22P020187/2022|EPI_ISL_19135640",
    "A/turkey/Italy/21VIR9869-4/2021|EPI_ISL_14760600",
    "A/Barnacle_goose/Scotland/058494/2021|EPI_ISL_13453561",
    "A/porpoise_/Sweden/SVA220712SZ0367/FB002184/O-2022|EPI_ISL_14810369",
    "A/goose/Czech_Republic/23458-4T/2021|EPI_ISL_9603775",
    "A/northern_gannet/France/22P020794/2022|EPI_ISL_19135815",
    "A/greylag_goose/Poland/MB503_21RS3290-18/2021|EPI_ISL_18245810",
    "A/Greylag_goose/Sweden/SVA220308SZ0382/FB000559/M-2022|EPI_ISL_11580255",
    "A/buzzard/Germany-SH/AI00547/2022|EPI_ISL_14570788",
    "A/goose/Italy/22VIR1520/2022|EPI_ISL_14761428",
    "A/Barnacle_goose/Scotland/003579/2022|EPI_ISL_13486700",
    "A/Great_black-backed_gull/Norway/2022-07-1235-4/2022|EPI_ISL_18901625",
    "A/White-tailed_eagle/Norway/2022-07-450-1/2022|EPI_ISL_18901680",
    "A/buzzard/Wales/055388/2021|EPI_ISL_13370343",
    "A/White-tailed_eagle/Norway/2022-07-984-1/2022|EPI_ISL_18901688",
    "A/barnacle_goose/Scotland/072143/2021|EPI_ISL_13370522",
    "A/turkey/Italy/21VIR9996/2021|EPI_ISL_14760615",
    "A/Canada_goose/Luxembourg/21237199/2021|EPI_ISL_13026105",
    "A/peacock/Italy/22VIR1953-2/2022|EPI_ISL_14761473",
    "A/northern_gannet/Sweden/SVA220630SZ0368/FB002052/N-2022|EPI_ISL_14985603",
    "A/chicken/Italy/22VIR1953-4/2022|EPI_ISL_14761496",
    "A/turkey/Italy/21VIR9906-2/2021|EPI_ISL_14760606",
    "A/duck/Poland/H1971_21RS3290-8/2021|EPI_ISL_18458062",
    "A/chicken/Italy/21VIR10146/2021|EPI_ISL_14760636",
    "A/turkey/England/075922/2021|EPI_ISL_13370620",
    "A/Turkey/Sweden/SVA211212SZ0001/FB301013-IP-2/M-2021|EPI_ISL_8338002",
    "A/chicken/Czech_Republic/22224-3T/2021|EPI_ISL_7224463",
    "A/layer/Ireland/034424_22VIR1325-21/2021|EPI_ISL_11259272",
    "A/laying_hen/Italy/21VIR11597-1/2021|EPI_ISL_14761065",
    "A/turkey/Israel/564/2021|EPI_ISL_12917757",
    "A/mute_swan/Northern_Ireland/17634/2021|EPI_ISL_13370521",
    "A/squared-tailed_nightjar/Nigeria/VRD22-305-08_23VIR8254-103/2022|EPI_ISL_18462531",
    "A/chicken/Czech_Republic/23404-4T/2021|EPI_ISL_7626513",
    "A/Common_Buzzard/Scotland/031887/2022|EPI_ISL_13969420",
    "A/chicken/Italy/22VIR1694-5/2022|EPI_ISL_14761451",
    "A/turkey/Italy/21VIR9609-2/2021|EPI_ISL_14760572",
    "A/Great_black-backed_gull/Norway/2022-07-2853-2/2022|EPI_ISL_18901669",
    "A/turkey/Italy/21VIR11587-1/2021|EPI_ISL_14761097",
    "A/Great_black-backed_gull/Norway/2022-07-2853-3/2022|EPI_ISL_18901670",
    "A/Gull/Norway/2022-07-3180-1/2022|EPI_ISL_18901673",
    "A/northern_gannet/Sweden/SVA220704SZ0158/FB002085/M-2022|EPI_ISL_14985874",
    "A/laying_hen/Italy/21VIR11169/2021|EPI_ISL_14760733",
    "A/greylag_goose/Northern_Ireland/17795/2021|EPI_ISL_13370556",
    "A/turkey/Italy/21VIR10988-1/2021|EPI_ISL_14760708",
    "A/Common_Buzzard/England/245285/2021|EPI_ISL_13370629",
    "A/great_crested_grebe/Netherlands/22001219-002/2022|EPI_ISL_10724258",
    "A/Goshawk/England/387180/2021|EPI_ISL_13370630",
    "A/Whooper_swan/Scotland/056219/2021|EPI_ISL_9029962",
    "A/turkey/Italy/21VIR9697-3_/2022|EPI_ISL_14760594",
    "A/northern_gannet/France/22P020191/2022|EPI_ISL_19135641",
    "A/barnacle_goose/Sweden/SVA211221SZ0510/FB004868/2021|EPI_ISL_9649405",
    "A/Chicken/Sweden/SVA211216SZ0001/FB306314/2021|EPI_ISL_9640022",
    "A/goose/Netherlands/22000310-001/2022|EPI_ISL_9616194",
    "A/greylag_goose/Netherlands/21039376-001/2021|EPI_ISL_7267254",
    "A/chicken/Italy/21VIR10969-3/2021|EPI_ISL_14760707",
    "A/laying_hen/Austria/22005859-001/2022|EPI_ISL_17785727",
    "A/common_buzzard/Netherlands/21039607-002/2021|EPI_ISL_7683101",
    "A/greylag_goose_/Sweden/SVA211118SZ0354/FB004497/I-2021|EPI_ISL_7055384",
    "A/swan/Poland/MB083_22VIR2515-8/2022|EPI_ISL_11922814",
    "A/northern_goshawk/Sweden/SVA220210SZ0305/FP000317/O-2022|EPI_ISL_10943657",
    "A/chicken/Scotland/130966/2022|EPI_ISL_15585888",
    "A/northern_gannet/France/22P019964/2022|EPI_ISL_19135633",
    "A/pheasant/Czech_Republic/25827-3/2021|EPI_ISL_9603939",
    "A/domestic_goose/England/008781/2022|EPI_ISL_13370698",
    "A/chicken/Italy/21VIR10522-13/2021|EPI_ISL_14760666",
    "A/chicken/Italy/IZSLT22VIR2562-1/2022|EPI_ISL_14761519",
    "A/great_black-backed_gull/Germany-SH/AI05236/2022|EPI_ISL_17693350",
    "A/Northern_goshawk/Norway/2023-07-77-1/2023|EPI_ISL_18963398",
    "A/turkey/Italy/21VIR9213-1/2021|EPI_ISL_14760555",
    "A/Greylag_goose/England/054503/2021|EPI_ISL_9029965",
    "A/chicken/Italy/21VIR10999-1/2021|EPI_ISL_14760722",
    "A/buzzard/Germany-TH/AI01791/2022|EPI_ISL_18007043",
    "A/turkey/Poland/H1977_21RS3290-10/2021|EPI_ISL_18458064",
    "A/northern_gannet/France/22P020063/2022|EPI_ISL_19135638",
    "A/Mute_swan/Wales/058560/2022|EPI_ISL_13969430",
    "A/swan/Austria/22046299-002/2022|EPI_ISL_17785728",
    "A/red_fox/Netherlands/22002993-003_/2022|EPI_ISL_10255587",
    "A/chicken/Italy/21VIR10388/2021|EPI_ISL_7733582",
    "A/chicken/Netherlands/21039108-001005/2021|EPI_ISL_6675739",
    "A/Eurasian_wigeon/Germany-SH/AI06143/2021|EPI_ISL_7753387",
    "A/barnacle_goose/Sweden/SVA211221SZ0510/FB004869/2021|EPI_ISL_9649612",
    "A/avian/Italy/21VIR11236/2021|EPI_ISL_14760742",
    "A/turkey/Italy/21VIR9372-1/2021|EPI_ISL_14760560",
    "A/turkey/Spain/140-38_22VIR2142-19/2022|EPI_ISL_11259291",
    "A/turkey/Denmark/24325-19/2021|EPI_ISL_19237886",
    "A/chicken/Italy/21VIR10723-1/2021|EPI_ISL_14760675",
    "A/Wild_goose/Italy/21VIR10193/2021|EPI_ISL_7733611",
    "A/swan/Austria/22030844-004/2022|EPI_ISL_17785733",
    "A/Grey_heron/Luxembourg/21243177/2021|EPI_ISL_13026236",
    "A/northern_gannet/France/22P020811/2022|EPI_ISL_19135826",
    "A/turkey/Italy/21VIR9769-4/2021|EPI_ISL_14760595",
    "A/Peregrine_falcon/Ireland/031358_21VIR9564-2/2021|EPI_ISL_18458054",
    "A/northern_gannet/France/22P019762/2022|EPI_ISL_19135619",
    "A/seagull/Italy/21VIR10481-2/2021|EPI_ISL_8882210",
    "A/turkey/Italy/21VIR8728-1/2021|EPI_ISL_14760542",
    "A/turkey/Italy/21VIR9211-1/2021|EPI_ISL_14760553",
    "A/turkey/Italy/21VIR11506/2022|EPI_ISL_14760902",
    "A/chicken/Netherlands/21037750-001005/2021|EPI_ISL_5942524",
    "A/duck/Poland/H1993_21RS3290-14/2021|EPI_ISL_18458069",
    "A/Northern_Gannet/Netherlands/5/2022|EPI_ISL_15267026",
    "A/laying_hen/Italy/22VIR203-3/2022|EPI_ISL_14761302",
    "A/domestic_duck/England/065343/2021|EPI_ISL_13370513",
    "A/northern_gannet/France/22P020314/2022|EPI_ISL_19135722",
    "A/turkey/Italy/21VIR11006-1/2021|EPI_ISL_14760724",
    "A/turkey/Italy/21VIR9949-9/2021|EPI_ISL_14760623",
    "A/turkey/Italy/21VIR9520-3/2021|EPI_ISL_17684608",
    "A/chicken/Italy/21VIR11072-1/2021|EPI_ISL_14760710",
    "A/turkey/Italy/21VIR9643-2/2021|EPI_ISL_14760579",
    "A/coot/England/388108/2022|EPI_ISL_13370894",
    "A/northern_gannet/Sweden/SVA220701SZ0245/FB002072/O-2022|EPI_ISL_14985765",
    "A/Goose/Scotland/036879/2022|EPI_ISL_13969422",
    "A/Great_skua/Iceland/2023AI04253/2022|EPI_ISL_18001846",
    "A/domestic_duck/Northern_Ireland/AI42_21VIR11350-5/2021|EPI_ISL_8810667",
    "A/nothern_gannet/France/22P019331/2022|EPI_ISL_14702899",
    "A/hawk/Italy/22VIR428-7/2022|EPI_ISL_14761349",
    "A/Greylag_goose/England/247696/2022|EPI_ISL_13969428",
    "A/Anser_anser/Belgium/00939-0011/2022|EPI_ISL_16676477",
    "A/mute_swan/Netherlands/22004560-001/2022|EPI_ISL_12474791",
    "A/Black-headed_gull/England/064433/2021|EPI_ISL_13370416",
    "A/barnacle_goose/Sweden/SVA220217SZ0388/FB00388/K-2022|EPI_ISL_10944835",
    "A/duck/Italy/22VIR1295-1/2022|EPI_ISL_14761417",
    "A/barnacle_goose/Denmark/24273-1.02/2021|EPI_ISL_18953100",
    "A/red_fox/Netherlands/22002636-002/2022|EPI_ISL_11019882",
    "A/Great_black-backed_gull/1/2022|EPI_ISL_9261742",
    "A/northern_gannet/France/22P019220/2022|EPI_ISL_19135604",
    "A/buzzard/Germany-MV/AI00055/2022|EPI_ISL_14570557",
    "A/turkey/Italy/21VIR11587-2/2021|EPI_ISL_14761107",
    "A/northern_goshawk/Netherlands/22000305-002/2022|EPI_ISL_9117228",
    "A/turkey/Italy/21VIR10097/2021|EPI_ISL_14760627",
    "A/goose/Czech_Republic/22608-1/2021|EPI_ISL_8515478",
    "A/goose/Czech_Republic/25322-179/2021|EPI_ISL_12150664",
    "A/turkey/Italy/21VIR9889-1/2021|EPI_ISL_14760613",
    "A/harbour_seal/Scotland/054948/2022|EPI_ISL_17465835",
    "A/turkey/Italy/21VIR9997/2021|EPI_ISL_14760618",
    "A/chicken/Czech_Republic/22224-4T/2021|EPI_ISL_7224472",
    "A/northern_gannet/France/22P019885/2022|EPI_ISL_19135630",
    "A/ibis/Egypt/RLQP-229S/2022|EPI_ISL_15614803",
    "A/turkey/Italy/21VIR8926-2/2021|EPI_ISL_14760546",
    "A/barnacle_goose/Germany-SH/AI06005/2021|EPI_ISL_7753251",
    "A/chicken/Denmark/24357-2/2021|EPI_ISL_18953102",
    "A/duck/Poland/H1982_21RS3290-12/2021|EPI_ISL_18458067",
    "A/swan/Romania/10656_22VIR2749-6/2022|EPI_ISL_11922820",
    "A/chicken/Poland/H1940-N/2021|EPI_ISL_6931288",
    "A/chicken/Tyumen/27-42V/2021|EPI_ISL_5463796",
    "A/Great_black-backed_gull/Norway/2022-07-1235-1/2022|EPI_ISL_18901622",
    "A/European_polecat/Netherlands/1/2022|EPI_ISL_13201074",
    "A/barnacle_goose/Scotland/064368/2021|EPI_ISL_13370417",
    "A/chicken/England/012967/2022|EPI_ISL_11561590",
    "A/turkey/Wales/065047/2021|EPI_ISL_13370511",
    "A/Gannet/Scotland/090910/2022|EPI_ISL_13969426",
    "A/barnacle_goose/England/294620/2021|EPI_ISL_13370561",
    "A/chicken/Czech_Republic/23404-2K/2021|EPI_ISL_7626484",
    "A/Northern_gannet/Norway/2022-07-1073-2/2022|EPI_ISL_18901610",
    "A/turkey/Italy/21VIR8816-1/2021|EPI_ISL_14760543",
    "A/domestic_duck/England/012973/2022|EPI_ISL_11561591",
    "A/chicken/Italy/21VIR10264/2021|EPI_ISL_14760645",
    "A/broiler/Italy/21VIR11886-1/2021|EPI_ISL_14761249",
    "A/chicken/Czech_Republic/61-1/2022|EPI_ISL_11327088",
    "A/turkey/Poland/H1944-N/2021|EPI_ISL_6935584",
    "A/chicken/Italy/21VIR10096/2021|EPI_ISL_14760625",
    "A/whooper_swan/England/055405/2021|EPI_ISL_13369741",
    "A/northern_gannet/Sweden/SVA220701SZ0245/FB002073/O-2022|EPI_ISL_14985873",
    "A/chicken/Czech_Republic/22224-2T/2021|EPI_ISL_7224437",
    "A/greylag_goose/Netherlands/21037497-001/2021|EPI_ISL_6328001",
    "A/turkey/Italy/21VIR10737-1/2021|EPI_ISL_14760690",
    "A/chicken/Italy/21VIR10252/2021|EPI_ISL_14760635",
    "A/chicken/Netherlands/21038165-006010/2021|EPI_ISL_13360259",
    "A/chicken/Italy/21VIR9133-20/2021|EPI_ISL_7733617",
    "A/mute_swan/Czech_Republic/2755/2022|EPI_ISL_12150661",
    "A/turkey/Poland/H1913-T1/2021|EPI_ISL_6930238",
    "A/chicken/Northern_Ireland/AI48_22VIR584-11/2021|EPI_ISL_18458591",
    "A/goose/Czech_Republic/25322-229/2021|EPI_ISL_12325995",
    "A/chicken/Italy/21VIR9944-2/2021|EPI_ISL_14760622",
    "A/buzzard/Scotland/043118/2022|EPI_ISL_13453565",
    "A/northern_goshawk/Germany-BB/AI01628/2022|EPI_ISL_14573032",
    "A/northern_gannet/Sweden/SVA220701SZ0245/FB002065/N-2022|EPI_ISL_14985606",
    "A/northern_gannet/Germany-SH/AI03498/2022|EPI_ISL_14573031",
    "A/laying_hen/Italy/21VIR11595-5/2021|EPI_ISL_14761052",
    "A/Harbour_Seal/Scotland/162925/2022|EPI_ISL_16904554",
    "A/chicken/Italy/21VIR9212-1/2021|EPI_ISL_14760552",
    "A/northern_gannet/Germany-NI/AI04234/2022|EPI_ISL_16122154",
    "A/laying_hen/Italy/21VIR11502/2021|EPI_ISL_8882173",
    "A/chicken/Italy/21VIR10469/2021|EPI_ISL_14760659",
    "A/pheasant/Italy/21VIR10359/2021|EPI_ISL_14760637",
    "A/northern_gannet/France/22P020485/2022|EPI_ISL_19135763",
    "A/turkey/Italy/21VIR9605-4/2021|EPI_ISL_14760573",
    "A/Tyto_alba/Belgium/334_0012/2021|EPI_ISL_14388346",
    "A/greylag_goose/Netherlands/21038417-001/2021|EPI_ISL_7952126",
    "A/turkey/Italy/21VIR11495-5/2021|EPI_ISL_14760967",
    "A/grey_heron/Croatia/132-a/2021|EPI_ISL_18361501",
    "A/chicken/Scotland/064936/2021|EPI_ISL_13370510",
    "A/mute_swan/Netherlands/21039291-002/2021|EPI_ISL_7996369",
    "A/domestic_duck/Germany-SH/AI06221/2021|EPI_ISL_14571443",
    "A/chicken/Italy/21VIR9507-7/2021|EPI_ISL_14760567",
    "A/chicken/Italy/21VIR10463-1/2021|EPI_ISL_14760665",
    "A/Nothern_gannet/France/22P019694/2022|EPI_ISL_18445751",
    "A/turkey/Italy/21VIR10095/2021|EPI_ISL_14760628",
    "A/barnacle_goose/Netherlands/22000551-002/2022|EPI_ISL_9304992",
    "A/greylag_goose/Netherlands/21038565-002/2021|EPI_ISL_6761021",
    "A/Guinea_fowl/Italy/21VIR10216/2021|EPI_ISL_14760632",
    "A/Common_buzzard/Scotland/073099/2021|EPI_ISL_13486819",
    "A/common_buzzard/Denmark/24271-1.02/2021-10-25|EPI_ISL_11798571",
    "A/wild_bird/Netherlands/22002833-002/2022|EPI_ISL_10255580",
    "A/turkey/Germany-MV/AI06035/2021|EPI_ISL_7753290",
    "A/turkey/Italy/21VIR11507/2021|EPI_ISL_12028893",
    "A/chicken/Czech_Republic/23589-3/2021|EPI_ISL_7626532",
    "A/northern_gannet/Spain/22VIR12774-17_3389-5/2022|EPI_ISL_16507100",
    "A/gannet/Ireland/025964_22VIR12757-2/2022|EPI_ISL_20077606",
    "A/goose/Chelyabinsk/1341-3/2021|EPI_ISL_16838578",
    "A/swan/Italy/22VIR1560/2022|EPI_ISL_14761383",
    "A/northern_gannet/France/22P020770/2022|EPI_ISL_19135796",
    "A/quail/Italy/21VIR10356/2021|EPI_ISL_14760639",
    "A/turkey/Italy/21VIR10001/2021|EPI_ISL_14760619",
    "A/chicken/Netherlands/21038669-001005/2021|EPI_ISL_6761017",
    "A/buzzard/Ireland/000656_22VIR1325-12/2022|EPI_ISL_11259263",
    "A/duck/Chelyabinsk/1379-1/2021|EPI_ISL_16838577",
    "A_common_buzzard_Denmark_24962-1.01_2021|EPI_ISL_18788584") # keep main clade

seq_C2 <- seq_C[(str_extract(names(seq_C), "EPI_[A-Z]+_\\d+") %in% str_extract(keep_C2, "EPI_[A-Z]+_\\d+"))] # extract through EPI IDs
rm_C2_clade <- names(seq_C[!(str_extract(names(seq_C), "EPI_[A-Z]+_\\d+") %in% str_extract(keep_C2, "EPI_[A-Z]+_\\d+"))])

## Remove: No exact collection date
meta_BE$sampling_date <- dmy(meta_BE$sampling_date)
meta_gisaid$collection_date <- dmy(meta_gisaid$collection_date)

rm_date_id <- meta_gisaid[which(is.na(meta_gisaid$collection_date)), "id_gisaid"]

rm_A_date <- names(seq_A)[str_extract(names(seq_A), "EPI_[A-Z]+_\\d+") %in% rm_date_id]
rm_C2_date <- names(seq_C2)[str_extract(names(seq_C2), "EPI_[A-Z]+_\\d+") %in% rm_date_id]
rm_AB_date <- names(seq_AB)[str_extract(names(seq_AB), "EPI_[A-Z]+_\\d+") %in% rm_date_id]
rm_BB_date <- names(seq_BB)[str_extract(names(seq_BB), "EPI_[A-Z]+_\\d+") %in% rm_date_id]

seq_A <- seq_A[!(str_extract(names(seq_A), "EPI_[A-Z]+_\\d+") %in% rm_date_id)]
seq_C2 <- seq_C2[!(str_extract(names(seq_C2), "EPI_[A-Z]+_\\d+") %in% rm_date_id)]
seq_AB <- seq_AB[!(str_extract(names(seq_AB), "EPI_[A-Z]+_\\d+") %in% rm_date_id)]
seq_BB <- seq_BB[!(str_extract(names(seq_BB), "EPI_[A-Z]+_\\d+") %in% rm_date_id)]

## Remove: Reassortment 
rm_C2_re <- c("A/domestic_duck/England/074574/2021|EPI_ISL_13370616",
    "A/mute_swan/Poland/MB008-22_22VIR5675-16/2022|EPI_ISL_18245597", 
    "A/European_polecat/Netherlands/1/2022|EPI_ISL_13201074",
    "A/greylag_goose/Poland/MB503_21RS3290-18/2021|EPI_ISL_18245810",
    "A/Northern_gannet/Norway/2022-07-1370-1/2022|EPI_ISL_18901633") # no recombination signal found

rm_AB_re <- c("A/hering_gull/Austria/23006560-001/2023|EPI_ISL_18058693",
  "A/canada_goose/England/119462/2022|EPI_ISL_16384170",
  "A/Domestic_Duck/England/137313/2023|EPI_ISL_18668521",
  "A/Domestic_Duck/England/137332/2023|EPI_ISL_18668522",
  "A/Mallard/Netherlands/19/2022|EPI_ISL_15581789",
  "A/Great_black-backed_Gull/Netherlands/3/2022|EPI_ISL_11560325",
  "A/domestic_duck/England/133434/2022|EPI_ISL_15586129",
  "A/chicken/Czech_Republic/22968_orig/2022|EPI_ISL_16613718",
  "A/Harris_Hawk/Scotland/014071/2023|EPI_ISL_16956028")

rm_BB_re <- c("A/black-headed_gull/Austria/23040809-002/2023|EPI_ISL_17716066",
  "A/European_Herring_Gull/Netherlands/9/2022|EPI_ISL_14163714",
  "A/black-headed_gull/Austria/23031105-002/2023|EPI_ISL_17514165",
  "A/peregrine_falcon/Germany-NI/2023AI07594/2023|EPI_ISL_18432789",
  "A/black-headed_gull/Austria/23052993-001/2023|EPI_ISL_17716068",
  "A/black-headed_gull/Austria/23045023-001/2023|EPI_ISL_17716067",
  "A/black-headed_gull/Austria/23024615/2023|EPI_ISL_17514176",
  "A/black-headed_gull/Austria/23029679/2023|EPI_ISL_17514133",
  "A/black-headed_gull/Austria/23056483-005/2023|EPI_ISL_17716072",
  "A/black-headed_gull/Austria/23025916-001/2023|EPI_ISL_17716063",
  "A/black-headed_gull/Austria/23033899-002/2023|EPI_ISL_17716065",
  "A/goshawk/Austria/23048718-001/2023|EPI_ISL_17716088",
  "A/black-headed_gull/Austria/23037203-001/2023|EPI_ISL_17514132",
  "A/black-headed_gull/Austria/23030829-001/2023|EPI_ISL_17716064",
  "A/black-headed_gull/Austria/23060170-008/2023|EPI_ISL_17716077",
  "A/chicken/England/029912/2023|EPI_ISL_17465831")
  
seq_C2 <- seq_C2[!(str_extract(names(seq_C2), "EPI_[A-Z]+_\\d+") %in% str_extract(rm_C2_re, "EPI_[A-Z]+_\\d+"))]
seq_AB <- seq_AB[!(str_extract(names(seq_AB), "EPI_[A-Z]+_\\d+") %in% str_extract(rm_AB_re, "EPI_[A-Z]+_\\d+"))]
seq_BB <- seq_BB[!(str_extract(names(seq_BB), "EPI_[A-Z]+_\\d+") %in% str_extract(rm_BB_re, "EPI_[A-Z]+_\\d+"))]

## Remove: Unsupported ambiguity code by SplitsTree (B)
rm_BB_code <- "A/Domestic_Duck/England/062587/2023|EPI_ISL_17710026"
seq_BB <- seq_BB[!(str_extract(names(seq_BB), "EPI_[A-Z]+_\\d+") %in% str_extract(rm_BB_code, "EPI_[A-Z]+_\\d+"))]

## Remove: IDs
rm_A <- c()
rm_C <- c(rm_C2_clade, rm_C2_re)
rm_AB <- c(rm_AB_re)
rm_BB <- c(rm_BB_re, rm_BB_code)

rm_all <- c(rm_A, rm_C, rm_AB, rm_BB, rm_date_id)

## Write sequence files
'write.fasta(sequences = seq_A, names = names(seq_A), file.out = here(base_dir, "data", "sequences", "new", "EA-2020-A_final.fasta"))
write.fasta(sequences = seq_C2, names = names(seq_C2), file.out = here(base_dir, "data", "sequences", "new", "EA-2020-C2_final.fasta"))
write.fasta(sequences = seq_AB, names = names(seq_AB), file.out = here(base_dir, "data", "sequences", "new", "EA-2021-AB_final.fasta"))
write.fasta(sequences = seq_BB, names = names(seq_BB), file.out = here(base_dir, "data", "sequences", "new", "EA-2022-BB_final.fasta"))'

## Filter metadata
meta_BE <- meta_BE[!(meta_BE$id_gisaid %in% str_extract(rm_all, "EPI_[A-Z]+_\\d+")),] # none removed
meta_gisaid <- meta_gisaid[!(meta_gisaid$id_gisaid %in% str_extract(rm_all, "EPI_[A-Z]+_\\d+")),]

'write.csv(meta_BE, here(base_dir, "data", "metadata", "meta_BE_final.csv"), row.names = FALSE, quote = FALSE)
write.csv(meta_gisaid, here(base_dir, "data", "metadata", "meta_gisaid_final.csv"), row.names = FALSE, quote = FALSE)'

meta_BE_miss <- meta_BE[which(!meta_BE$id_gisaid %in% meta_gisaid$id_gisaid), c("id_gisaid", "id_name", "sampling_date", "location_province")]
colnames(meta_BE_miss) <- c("id_gisaid", "id_name", "collection_date", "location_province")
meta_BE_miss$location_continent  <- "Europe"
meta_BE_miss$location_country <- "Belgium"

meta_all <- rbind(meta_gisaid, meta_BE_miss)

#write.csv(meta_all, here(base_dir, "data", "metadata", "meta_all.csv"), row.names = FALSE)

meta_filter_func <- function(seqs, meta_df) {
    seq_names <- names(seqs)
    seq_ids <- str_extract(seq_names, "EPI_[A-Z]+_\\d+")
    lookup <- tibble(id_gisaid = seq_ids, seq_name  = seq_names)
    meta_df %>% filter(id_gisaid %in% seq_ids) %>% left_join(lookup, by = "id_gisaid") %>% relocate(seq_name, .before = 1)
}

meta_A <- meta_filter_func(seq_A, meta_all)
meta_C <- meta_filter_func(seq_C, meta_all)
meta_AB <- meta_filter_func(seq_AB, meta_all)
meta_BB <- meta_filter_func(seq_BB, meta_all)

'write.csv(meta_A, here(base_dir, "data", "metadata", "meta_EA-2020-A.csv"), row.names = FALSE, quote = FALSE)
write.csv(meta_C, here(base_dir, "data", "metadata", "meta_EA-2020-C.csv"), row.names = FALSE, quote = FALSE)
write.csv(meta_AB, here(base_dir, "data", "metadata", "meta_EA-2021-AB.csv"), row.names = FALSE, quote = FALSE)
write.csv(meta_BB, here(base_dir, "data", "metadata", "meta_EA-2022-BB.csv"), row.names = FALSE, quote = FALSE)'

# Epidemiological data: Occurrence data from internal Sciensano NeoLIMS database

epi_data_raw <- read_csv(here(base_dir, "data", "metadata", "occurrences_public.csv"))

epi_data <- epi_data_raw %>% 
    filter(Result_AFLoccurence == "POSITIVE") %>% 
    filter(Result_H5H7Other == "H5") %>% 
    mutate(species_group = case_when(
        animal_type == "Wild" ~ "wild_birds",
        animal_type == "Farm" ~ "poultry",
        animal_type == "Hobby" ~ "poultry",
        animal_type == "Hobby?" ~ "poultry",
        animal_type == "NoBird" ~ "no_bird",
        animal_type == "Zoo" ~ "captive_birds",
        TRUE ~ "other"
    ))

epi_data <- epi_data %>%
    mutate(decimal_date = lubridate::decimal_date(lubridate::dmy(final_date))) %>% 
    mutate(species_group = factor(species_group, levels = c("wild_birds", "poultry", "captive_birds", "no_bird")))

table(epi_data$species_group)

# TreeTime

## Date input
treetime_A <- meta_A %>% dplyr::select(seq_name, collection_date) %>% rename(name = seq_name, date = collection_date) %>% mutate(name = str_replace_all(name, "[()]", "_"))
treetime_C <- meta_C %>% dplyr::select(seq_name, collection_date) %>% rename(name = seq_name, date = collection_date) %>% mutate(name = str_replace_all(name, "[()]", "_"))
treetime_AB <- meta_AB %>% dplyr::select(seq_name, collection_date) %>% rename(name = seq_name, date = collection_date) %>% mutate(name = str_replace_all(name, "[()]", "_"))
treetime_BB <- meta_BB %>% dplyr::select(seq_name, collection_date) %>% rename(name = seq_name, date = collection_date) %>% mutate(name = str_replace_all(name, "[()]", "_"))

'write.csv(treetime_A, here(base_dir, "data", "metadata", "treetime_EA-2020-A.csv"), row.names = FALSE, quote = FALSE)
write.csv(treetime_C, here(base_dir, "data", "metadata", "treetime_EA-2020-C2.csv"), row.names = FALSE, quote = FALSE)
write.csv(treetime_AB, here(base_dir, "data", "metadata", "treetime_EA-2021-AB.csv"), row.names = FALSE, quote = FALSE)
write.csv(treetime_BB, here(base_dir, "data", "metadata", "treetime_EA-2022-BB.csv"), row.names = FALSE, quote = FALSE)'

## Run TreeTime
# Use system() to run directly from code instead of from terminal
# treetime --time-marginal false --outdir final --tree EA-2022-BB.tre --aln EA-2022-BB_alignment.fasta --dates treetime_EA-2022-BB.csv

# DTA in BEAST

## Metadata
DTA_A <- meta_A %>% dplyr::select(seq_name, collection_date, location_country) %>% rename(id = seq_name, date = collection_date, country = location_country) %>% 
    mutate(location = if_else(country == "Belgium", "Belgium", "other")) %>% mutate(id = str_replace_all(id, "[()]", "_"))
DTA_C <- meta_C %>% dplyr::select(seq_name, collection_date, location_country) %>% rename(id = seq_name, date = collection_date, country = location_country) %>% 
    mutate(location = if_else(country == "Belgium", "Belgium", "other")) %>% mutate(id = str_replace_all(id, "[()]", "_"))
DTA_AB <- meta_AB %>% dplyr::select(seq_name, collection_date, location_country) %>% rename(id = seq_name, date = collection_date, country = location_country) %>% 
    mutate(location = if_else(country == "Belgium", "Belgium", "other")) %>% mutate(id = str_replace_all(id, "[()]", "_")) 
DTA_BB <- meta_BB %>% dplyr::select(seq_name, collection_date, location_country) %>% rename(id = seq_name, date = collection_date, country = location_country) %>% 
    mutate(location = if_else(country == "Belgium", "Belgium", "other")) %>% mutate(id = str_replace_all(id, "[()]", "_"))

'write.table(DTA_A, here(base_dir, "data", "metadata", "BEAST_EA-2020-A.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(DTA_C, here(base_dir, "data", "metadata", "BEAST_EA-2020-C2.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(DTA_AB, here(base_dir, "data", "metadata", "BEAST_EA-2021-AB.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(DTA_BB, here(base_dir, "data", "metadata", "BEAST_EA-2022-BB.txt"), row.names = FALSE, quote = FALSE, sep = "\t")'

## Write XML from template
genotype_vec <- c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB") 

chain_length <- "200000"
DTA_template_base <- scan(here(base_dir, "data", "BEAST", "BEAST_DTA_template.xml"), what = "", sep = "\n", quiet = TRUE, blank.lines.skip = FALSE)
DTA_template_base <- gsub('chainLength="[0-9]+"', paste0('chainLength="', chain_length, '"'), DTA_template_base)

treetime_data <- list()

for (i in genotype_vec) {

    # IMPORTANT: fresh template per genotype
    DTA_template <- gsub("TEMPLATE", paste0(i, "_DTA"), DTA_template_base)

    treetime_tree <- read.nexus(here(base_dir, "data", "treetime", i, "final", "timetree.nexus"))
    treetime_tree <- write.tree(treetime_tree)
    DTA_meta <- read.table(here(base_dir, "data", "BEAST", paste0("BEAST_", i, ".txt")), sep = "\t", header = TRUE)

    treetime_data[[i]] <- list(
        treetime_tree = treetime_tree,
        treetime_meta = read.csv(here(base_dir, "data", "treetime", i, paste0("treetime_", i, ".csv"))),
        DTA_meta = DTA_meta
    )

    out_dir <- here(base_dir, "data", "BEAST", i, "DTA_run")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    sink(file = here(out_dir, paste0("BEAST_DTA_", i, ".xml")))

    for (j in seq_along(DTA_template)) {
        cat(DTA_template[j], "\n")

        if (grepl("<taxa id=\"taxa\">", DTA_template[j])) {
            for (k in seq_len(nrow(DTA_meta))) {
                cat("\t\t<taxon id=\"", DTA_meta[k, "id"], "\">\n", sep = "")
                cat("\t\t\t<date value=\"", decimal_date(ymd(as.character(DTA_meta[k, "date"]))), "\" direction=\"forwards\" units=\"years\"/>\n", sep = "")
                cat("\t\t\t<attr name=\"location\">\n", sep = "")
                cat("\t\t\t\t", DTA_meta[k, "location"], "\n", sep = "")
                cat("\t\t\t</attr>\n", sep = "")
                cat("\t\t</taxon>\n", sep = "")
            }
        }

        if (grepl("<alignment id=\"alignment\" dataType=\"nucleotide\">", DTA_template[j])) {
            for (k in seq_len(nrow(DTA_meta))) {
                cat("\t\t<sequence>\n", sep = "")
                cat("\t\t\t<taxon idref=\"", DTA_meta[k, "id"], "\"/>\n", sep = "")
                cat("\t\t\t\tNNNN\n", sep = "")
                cat("\t\t</sequence>\n", sep = "")
            }
        }

        if (grepl("<newick id=\"startingTree\">", DTA_template[j])) {
            cat("\t\t", treetime_tree, "\n", sep = "")
        }
    }

    sink(NULL)
}

## Run BEAST
for(l in genotype_vec) {
    for(prefix in c("1_", "2_")) { # each XML-file is run twice
        local({
            wd <- getwd()
            on.exit(setwd(wd), add = TRUE)

            run_dir <- here(base_dir, "data", "BEAST", l, "DTA_run")
            xml <- paste0("BEAST_DTA_", l, ".xml")

            setwd(run_dir)

            system2(
                command = "beast",
                args = c("-overwrite", "-prefix", prefix ,xml), 
                stdout = file.path(run_dir, paste0(prefix, "beast_stdout", ".log")),
                stderr = file.path(run_dir, paste0(prefix, "beast_stderr", ".log")),
                wait = TRUE
            )
        })
    }
}

## Combine log-files

genotype_vec

for (x in genotype_vec) {
  local({
    wd <- getwd()
    on.exit(setwd(wd), add = TRUE)

    run_dir <- here(base_dir, "data", "BEAST", x, "DTA_run")
    setwd(run_dir)

    ## 1) Combine LOG files
    ## logcombiner -burnin 20000 1_x_DTA1.log 2_x_DTA1.log x_DTA.log
    system2(
      command = "logcombiner",
      args = c(
        "-burnin", "20000",
        paste0("1_", x, "_DTA1.log"),
        paste0("2_", x, "_DTA1.log"),
        paste0(x, "_DTA.log")
      ),
      stdout = file.path(run_dir, paste0(x, "_logcombiner_log_stdout.log")),
      stderr = file.path(run_dir, paste0(x, "_logcombiner_log_stderr.log")),
      wait = TRUE
    )

    ## 2) Combine TREE files
    ## logcombiner -trees -burnin 20000 1_x_DTA.trees 2_x_DTA.trees x_DTA.trees
    system2(
      command = "logcombiner",
      args = c(
        "-trees",
        "-burnin", "20000",
        paste0("1_", x, "_DTA.trees"),
        paste0("2_", x, "_DTA.trees"),
        paste0(x, "_DTA.trees")
      ),
      stdout = file.path(run_dir, paste0(x, "_logcombiner_trees_stdout.log")),
      stderr = file.path(run_dir, paste0(x, "_logcombiner_trees_stderr.log")),
      wait = TRUE
    )

    ## 3) Annotate trees using HIPSTR
    ## treeannotator -burnin 0 x_DTA.trees x_DTA.tree
    system2(
      command = "treeannotator",
      args = c(
        "-burnin", "0",
        paste0(x, "_DTA.trees")
        paste0(x, "_DTA.tree")
      ),
      stdout = file.path(run_dir, paste0(x, "_treeannotator_stdout.log")),
      stderr = file.path(run_dir, paste0(x, "_treeannotator_stderr.log")),
      wait = TRUE
    )
  })
}

# Estimating number of distinct introduction events
## Note: Only correct calculation if there are no posterior trees with 0 introductions

intro_summary_df <- tibble(
    genotype = character(),

    intro_median = numeric(),
    intro_95HPD_lower = numeric(),
    intro_95HPD_upper = numeric(),

    cluster_n1 = numeric(),
    cluster_n1_95HPD_lower = numeric(),
    cluster_n1_95HPD_upper = numeric(),

    clade_size_median = numeric(),
    clade_size_mean_95HPD_lower = numeric(),
    clade_size_mean_95HPD_upper = numeric(),

    clade_size_max = numeric(),
    clade_size_max_95HPD_lower = numeric(),
    clade_size_max_95HPD_upper = numeric()
)

genotype_vec

for(genotype_index in seq_along(genotype_vec)) {

    genotype <- genotype_vec[genotype_index]

    cat("Reading trees for", genotype, "...\n")
    timer <- system.time({
        trees <- read.beast(here(base_dir, "data", "BEAST", genotype, "DTA_run", paste0(genotype, "_DTA.trees")))
    })
    cat("Finished reading trees for", genotype, "in", timer["elapsed"], "seconds.\n")

    intro_counts <- numeric(length(trees))
    cluster_n1_abs_vec <- numeric(length(trees))
    cluster_n1_vec <- numeric(length(trees))
    mean_clade_size_vec <- numeric(length(trees))
    max_clade_size_vec <- numeric(length(trees))

    for (tree_index in seq_along(trees)) {
        tree <- trees[[tree_index]]@phylo
        tree_annotation <- as.data.frame(trees[[tree_index]]@data)

        tree_edge <- tree$edge
        tree_parent_node <- tree_edge[,1]
        tree_child_node <- tree_edge[,2]

        node_loc <- setNames(as.character(tree_annotation$location), as.character(tree_annotation$node))
        start_location <- node_loc[as.character(tree_parent_node)]
        end_location <- node_loc[as.character(tree_child_node)]

        intro_edges <- which(end_location == "Belgium" & start_location == "other")
        intro_counts[tree_index] <- length(intro_edges)

        tips_count <- length(tree$tip.label)
        cluster_n1 <- sum(tree_child_node[intro_edges] <= tips_count)
        cluster_n1_abs_vec[tree_index] <- cluster_n1
        cluster_n1_vec[tree_index] <- cluster_n1 / length(intro_edges)

        clade_sizes <- numeric(length(intro_edges))
        for(intro_index in seq_along(intro_edges)) {
            child_node <- tree_child_node[intro_edges[intro_index]]

            if (child_node <= tips_count) { # introduction leads to tip
                clade_sizes[intro_index] <- 1 
            } else {
                tip_labels <- extract.clade(tree, child_node)$tip.label
                tip_idx <- match(tip_labels, tree$tip.label)
                clade_sizes[intro_index] <- sum(node_loc[as.character(tip_idx)] == "Belgium")
            }
        }
        mean_clade_size_vec[tree_index] <- mean(clade_sizes)
        max_clade_size_vec[tree_index] <- max(clade_sizes)
    }

    intro_hpd <- hdi(intro_counts)[1:2]
    cluster_n1_hpd <- hdi(cluster_n1_vec)[1:2]
    cluster_n1_abs_hpd <- hdi(cluster_n1_abs_vec)[1:2]
    mean_clade_hpd <- hdi(mean_clade_size_vec)[1:2]
    max_clade_hpd <- hdi(max_clade_size_vec)[1:2]

    intro_summary_df <- bind_rows(
        intro_summary_df,
        tibble(
            genotype = genotype,

            intro_median = median(intro_counts),
            intro_95HPD_lower = intro_hpd[1],
            intro_95HPD_upper = intro_hpd[2],

            cluster_n1_abs = median(cluster_n1_abs_vec),
            cluster_n1_abs_95HPD_lower = cluster_n1_abs_hpd[1],
            cluster_n1_abs_95HPD_upper = cluster_n1_abs_hpd[2],

            cluster_n1 = median(cluster_n1_vec),
            cluster_n1_95HPD_lower = cluster_n1_hpd[1],
            cluster_n1_95HPD_upper = cluster_n1_hpd[2],

            clade_size_median = median(mean_clade_size_vec),
            clade_size_mean_95HPD_lower = mean_clade_hpd[1],
            clade_size_mean_95HPD_upper = mean_clade_hpd[2],
            
            clade_size_max = median(max_clade_size_vec),
            clade_size_max_95HPD_lower = max_clade_hpd[1],
            clade_size_max_95HPD_upper = max_clade_hpd[2]
        )
    )
    cat(genotype, "finished.\n")
}
intro_summary_df
write.csv(intro_summary_df, here(base_dir, "data", "metadata", "introduction_count.csv"), row.names = FALSE, quote = FALSE)

intro_summary_df <- read.csv(here(base_dir, "data", "metadata", "introduction_count.csv"), sep = ",", header = TRUE)
    
# Define clades for continous phylogeography 
## clusters1: list with vectors of tip labels for each introduced clade
## clusters2: clusters2[[genotype]][[cluster_index]]

clusters2 <- list()
intro_BE_list <- list()

for (genotype_index in seq_along(genotype_vec)) {

    genotype <- genotype_vec[genotype_index]

    tree_read <- read.beast(here(base_dir, "data", "BEAST", genotype, "DTA_run", paste0(genotype, "_DTA.tree")))
    tree <- tree_read@phylo
    tree_annotation <- as.data.frame(tree_read@data)

    node_loc <- setNames(as.character(tree_annotation$location.rate), as.character(tree_annotation$node))    

    tree_edge <- tree$edge
    tree_parent_node <- tree_edge[, 1]
    tree_child_node <- tree_edge[, 2]

    start_location <- node_loc[as.character(tree_parent_node)]
    end_location <- node_loc[as.character(tree_child_node)]
    intro_edges <- which(end_location == "Belgium" & start_location == "other")

    intro_BE_list[[genotype]] <- intro_edges

    tips_count <- length(tree$tip.label)
    clusters1 <- vector("list", length(intro_edges)) 

    for (intro_index in seq_along(intro_edges)) {
        child_node <- tree_child_node[intro_edges[intro_index]]
        clusters1[[intro_index]] <- if (child_node <= tips_count) tree$tip.label[child_node] else extract.clade(tree, node = child_node)$tip.label
    }

    clusters2[[genotype]] <- vector("list", length(clusters1))

    for (cluster_index in seq_along(clusters1)) {
        tip_id <- clusters1[[cluster_index]]
        if (length(tip_id) == 0) {
            clusters2[[genotype]][[cluster_index]] <- data.frame(tip_id = character(), sampling_date = character(), latitude = numeric(), longitude = numeric(), species_group = character(), species_host = character(), stringsAsFactors = FALSE)
        } else {
            tip_gisaid <- str_extract(tip_id, "EPI_[A-Z]+_\\d+")
            tab <- tibble(tip_id = tip_id, id_gisaid = tip_gisaid) %>%
                left_join(meta_BE, by = "id_gisaid") %>%
                dplyr::select(tip_id, sampling_date, latitude, longitude, species_group, species_host)
            clusters2[[genotype]][[cluster_index]] <- as.data.frame(tab, stringsAsFactors = FALSE)
        }
    }
}

#saveRDS(clusters2, here(base_dir, "data", "metadata", "intro_clusters.rds"))

for (genotype in names(clusters2)) { # decimal date
    for (i in seq_along(clusters2[[genotype]])) {
        clusters2[[genotype]][[i]]$sampling_date <- lubridate::decimal_date(ymd(clusters2[[genotype]][[i]]$sampling_date))
    }
}

sapply(clusters2, length)

# RRW in BEAST

## Write XML from template
genotype_vec
RRW_template_base <- scan(here(base_dir, "data", "BEAST", "BEAST_RRW_template.xml"), what = "", sep = "\n", quiet = TRUE, blank.lines.skip = FALSE)

chain_length <- "1000000000"
RRW_template_base <- gsub('chainLength="[0-9]+"', paste0('chainLength="', chain_length, '"'), RRW_template_base)

for(genotype in genotype_vec) {
    RRW_template <- gsub("TEMPLATE", paste0(genotype, "_RRW"), RRW_template_base)
    clusters <- clusters2[[genotype]]

    if (genotype == "EA-2021-AB") { # subsampling by phylogenetic clustering
        clusters[[9]] <- subset(clusters[[9]], tip_id != "A/Branta_canadensis/Belgium/5177_0003/2022|EPI_ISL_14393671")
        clusters[[26]] <- subset(clusters[[26]], !tip_id %in% c("A/Gallus_gallus/Belgium/14313-001/2022|EPI_ISL_16547851", "A/Mustela_furo/Belgium/14312-0002/2022|EPI_ISL_16547666"))
    }

    intro_BE <- intro_BE_list[[genotype]]

    out_dir <- here(base_dir, "data", "BEAST", genotype, "RRW_run")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    tree_read <- read.beast(here(base_dir, "data", "BEAST", genotype, "DTA_run", paste0(genotype, "_DTA.tree")))
    tree <- tree_read@phylo

    sink(file = here(out_dir, paste0("BEAST_RRW_", genotype, ".xml")))

    for (i in 1:length(RRW_template)) {
        cat(RRW_template[i], "\n")

        if (grepl("Insert taxa blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {

                    tips_rm <- character(0)
                    if (genotype == "EA-2021-AB" & j == 9) {
                        tips_rm <- c(
                            "A/Branta_canadensis/Belgium/5177_0003/2022|EPI_ISL_14393671"
                        )
                    }
                    if (genotype == "EA-2021-AB" & j == 26) {
                        tips_rm <- c(
                            "A/Gallus_gallus/Belgium/14313-001/2022|EPI_ISL_16547851",
                            "A/Mustela_furo/Belgium/14312-0002/2022|EPI_ISL_16547666"
                        )
                    }

                    cat(paste0("\t<taxa id=\"taxa_", j, "\">", "\n"))
                    for (k in 1:dim(clusters[[j]])[1]) {
                        if (!is.na(clusters[[j]][k, "latitude"]) & !clusters[[j]][k, "tip_id"] %in% tips_rm) {
                            cat(paste0("\t\t<taxon id=\"", clusters[[j]][k, "tip_id"], "\">", "\n"))
                            cat(paste0("\t\t\t<date value=\"", clusters[[j]][k, "sampling_date"], "\" direction=\"forwards\" units=\"years\"/>", "\n"))
                            cat("\t\t\t<attr name=\"latitude\">\n")
                            cat(paste0("\t\t\t\t", clusters[[j]][k, "latitude"], "\n"))
                            cat("\t\t\t</attr>\n")
                            cat("\t\t\t<attr name=\"longitude\">\n")
                            cat(paste0("\t\t\t\t", clusters[[j]][k, "longitude"], "\n"))
                            cat("\t\t\t</attr>\n")
                            cat("\t\t\t<attr name=\"location\">\n")
                            cat(paste0("\t\t\t\t", clusters[[j]][k, "latitude"], " ", clusters[[j]][k, "longitude"], "\n"))
                            cat("\t\t\t</attr>\n")
                            cat("\t\t</taxon>\n")
                        }
                    }
                    cat("\t</taxa>", "\n")
                }
            }
        }

        if (grepl("Insert alignment blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
 
                    cat(paste0("\t<alignment id=\"alignment_", j, "\" dataType=\"nucleotide\">", "\n"))
                    for (k in 1:dim(clusters[[j]])[1]) {
                        if (!is.na(clusters[[j]][k, "latitude"]) & !clusters[[j]][k, "tip_id"] %in% tips_rm) {
                            cat("\t\t<sequence>\n")
                            cat(paste0("\t\t\t<taxon idref=\"", clusters[[j]][k, "tip_id"], "\"/>", "\n"))
                            cat("\t\t\tNNNN\n")
                            cat("\t\t</sequence>\n")
                        }
                    }
                    cat("\t</alignment>", "\n")
                }
            }
        }

        if (grepl("Insert pattern blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t<patterns id=\"patterns_", j, "\" from=\"1\" strip=\"false\">", "\n"))
                    cat(paste0("\t\t<alignment idref=\"alignment_", j, "\"/>", "\n"))
                    cat("\t</patterns>", "\n")
                }
            }
        }

        if (grepl("Insert starting tree blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {

                    tips_rm <- character(0)
                    if (genotype == "EA-2021-AB" & j == 9) {
                        tips_rm <- c(
                            "A/Branta_canadensis/Belgium/5177_0003/2022|EPI_ISL_14393671"
                        )
                    }
                    if (genotype == "EA-2021-AB" & j == 26) {
                        tips_rm <- c(
                            "A/Gallus_gallus/Belgium/14313-001/2022|EPI_ISL_16547851",
                            "A/Mustela_furo/Belgium/14312-0002/2022|EPI_ISL_16547666"
                        )
                    }

                    child_node <- tree$edge[intro_BE[j], 2]
                    tips_count <- length(tree$tip.label)
                    tre <- if (child_node <= tips_count) ape::keep.tip(tree, tree$tip.label[child_node]) else ape::extract.clade(tree, node = child_node)

                    tips <- clusters[[j]][!is.na(clusters[[j]][, "latitude"]) & !clusters[[j]][, "tip_id"] %in% tips_rm, "tip_id"]
                    tips_to_keep <- tre$tip.label[which(gsub("'", "", tre$tip.label) %in% tips)]
                    tre <- ape::keep.tip(tre, tips_to_keep)

                    write.tree(tre, here(out_dir, paste0("Clade_", j, ".tre")))
                    tre <- scan(here(out_dir, paste0("Clade_", j, ".tre")), what = "", sep = "\n", quiet = TRUE)
                    txt <- c("#NEXUS", "begin trees;", paste0("\ttree tree_1 = [&R] ", tre), "end;")
                    write(txt, here(out_dir, paste0("Clade_", j, ".tre")))

                    cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_", j, "\" fileName=\"Clade_", j, ".tre\">", "\n"))
                    cat(paste0("\t\t<taxa idref=\"taxa_", j, "\"/>", "\n"))
                    cat("\t</empiricalTreeDistributionModel>", "\n")
                }
            }
        }

        if (grepl("Insert tree model blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t<treeModel id=\"treeModel_", j, "\">", "\n"))
                    cat(paste0("\t\t<coalescentTree idref=\"startingTree_", j, "\"/>", "\n"))
                    cat("\t\t<rootHeight>", "\n")
                    cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_", j, "\"/>", "\n"))
                    cat("\t\t</rootHeight>", "\n")
                    cat("\t\t<nodeHeights internalNodes=\"true\">", "\n")
                    cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_", j, "\"/>", "\n"))
                    cat("\t\t</nodeHeights>", "\n")
                    cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">", "\n")
                    cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_", j, "\"/>", "\n"))
                    cat("\t\t</nodeHeights>", "\n")
                    cat("\t</treeModel>", "\n")
                }
            }
        }

        if (grepl("Insert arbitraryBranchRates blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t<arbitraryBranchRates id=\"location.diffusion.branchRates", j, "\">", "\n"))
                    cat(paste0("\t\t<treeModel idref=\"treeModel_", j, "\"/>", "\n"))
                    cat("\t\t<rates>", "\n")
                    cat(paste0("\t\t\t<parameter id=\"location.diffusion.rates", j, "\" lower=\"0.0\"/>", "\n"))
                    cat("\t\t</rates>", "\n")
                    cat("\t</arbitraryBranchRates>", "\n")
                }
            }
        }

        if (grepl("Insert distributionLikelihood blocks 1", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t<distributionLikelihood id=\"location.diffusion.prior", j, "\">", "\n"))
                    cat("\t\t<data>", "\n")
                    cat(paste0("\t\t\t<parameter idref=\"location.diffusion.rates", j, "\"/>", "\n"))
                    cat("\t\t</data>", "\n")
                    cat("\t\t<distribution>", "\n")
                    cat("\t\t\t<onePGammaDistributionModel>", "\n")
                    cat("\t\t\t\t<shape>", "\n")
                    cat("\t\t\t\t\t<parameter value=\"0.5\"/>", "\n")
                    cat("\t\t\t\t</shape>", "\n")
                    cat("\t\t\t</onePGammaDistributionModel>", "\n")
                    cat("\t\t</distribution>", "\n")
                    cat("\t</distributionLikelihood>", "\n")
                }
            }
        }

        if (grepl("Insert location.traitLikelihood blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t<multivariateTraitLikelihood id=\"location.traitLikelihood", j, "\" traitName=\"location\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">", "\n"))
                    cat("\t\t<multivariateDiffusionModel idref=\"location.diffusionModel\"/>", "\n")
                    cat(paste0("\t\t<treeModel idref=\"treeModel_", j, "\"/>"))
                    cat("\t\t<traitParameter>", "\n")
                    cat(paste0("\t\t\t<parameter id=\"leaf.location", j, "\"/>", "\n"))
                    cat("\t\t</traitParameter>", "\n")
                    cat("\t\t<conjugateRootPrior>", "\n")
                    cat("\t\t\t<meanParameter>", "\n")
                    cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>", "\n")
                    cat("\t\t\t</meanParameter>", "\n")
                    cat("\t\t\t<priorSampleSize>", "\n")
                    cat("\t\t\t\t<parameter value=\"0.000001\"/>", "\n")
                    cat("\t\t\t</priorSampleSize>", "\n")
                    cat("\t\t</conjugateRootPrior>", "\n")
                    cat(paste0("\t\t<arbitraryBranchRates idref=\"location.diffusion.branchRates", j, "\"/>", "\n"))
                    cat("\t</multivariateTraitLikelihood>", "\n")
                }
            }
        }

        if (grepl("Insert continuousDiffusionStatistic blocks 1", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t<continuousDiffusionStatistic id=\"location.diffusionRate", j, "\" greatCircleDistance=\"true\">", "\n"))
                    cat(paste0("\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood", j, "\"/>", "\n"))
                    cat("\t</continuousDiffusionStatistic>", "\n")
                }
            }
        }

        if (grepl("Insert scaleOperator blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">", "\n")
                    cat(paste0("\t\t\t<parameter idref=\"location.diffusion.rates", j, "\"/>", "\n"))
                    cat("\t\t</scaleOperator>", "\n")
                }
            }
        }

        if (grepl("Insert precisionGibbsOperator blocks", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat("\t\t<precisionGibbsOperator weight=\"2\">", "\n")
                    cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood", j, "\"/>", "\n"))
                    cat("\t\t\t<multivariateWishartPrior idref=\"location.precisionPrior\"/>", "\n")
                    cat("\t\t</precisionGibbsOperator>", "\n")
                }
            }
        }

        if (grepl("Insert distributionLikelihood blocks 2", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t\t\t\t<distributionLikelihood idref=\"location.diffusion.prior", j, "\"/>", "\n"))
                }
            }
        }

        if (grepl("Insert multivariateTraitLikelihood blocks 1", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood", j, "\"/>", "\n"))
                }
            }
        }

        if (grepl("Insert continuousDiffusionStatistic blocks 2", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"location.diffusionRate", j, "\"/>", "\n"))
                }
            }
        }

        if (grepl("Insert multivariateTraitLikelihood blocks 2", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood", j, "\"/>", "\n"))
                }
            }
        }

        if (grepl("<!-- Insert logTree blocks -->", RRW_template[i])) {
            for (j in 1:length(clusters)) {
                if ((dim(clusters[[j]])[1] >= 3) & (sum(!is.na(clusters[[j]][, "longitude"])) >= 3)) {
                    cat(paste0("\t\t<logTree id=\"treeFileLog", j, "\" logEvery=\"100000\" nexusFormat=\"true\" fileName=\"Clade_", j, ".trees\" sortTranslationTable=\"true\">", "\n"))
                    cat(paste0("\t\t\t<treeModel idref=\"treeModel_", j, "\"/>", "\n"))
                    cat("\t\t\t<joint idref=\"joint\"/>", "\n")
                    cat("\t\t\t<trait name=\"location\" tag=\"location\">", "\n")
                    cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood", j, "\"/>", "\n"))
                    cat("\t\t\t</trait>\n")
                    cat("\t\t\t<multivariateDiffusionModel idref=\"location.diffusionModel\"/>\n")
                    cat("\t\t\t<trait name=\"rate\" tag=\"location.rate\">\n")
                    cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"location.diffusion.branchRates", j, "\"/>\n"))
                    cat("\t\t\t</trait>\n")
                    cat("\t\t</logTree>\n")
                }
            }
        }
    }
    sink(NULL)
}

cluster_id_list <- split(cluster_size_table$cluster_id, cluster_size_table$genotype)

clusters_rrw <- lapply(clusters2, function(clusters) { 
    clusters <- lapply(seq_along(clusters), function(cluster_id) {
        cl <- clusters[[cluster_id]] %>%
            mutate(cluster_id = cluster_id, .before = 1)
        cl <- cl[!is.na(cl$latitude), ]
        cl <- cl[!cl$tip_id %in% c(
            "A/Branta_canadensis/Belgium/5177_0003/2022|EPI_ISL_14393671",
            "A/Gallus_gallus/Belgium/14313-001/2022|EPI_ISL_16547851",
            "A/Mustela_furo/Belgium/14312-0002/2022|EPI_ISL_16547666"
        ), ]
        cl
    })
    clusters[sapply(clusters, nrow) >= 3]
})

## EA-2021-AB: Adding jitter for cluster 24
xml_AB_rrw <- scan(here(base_dir, "data", "BEAST", "EA-2021-AB", "RRW_run", "BEAST_RRW_EA-2021-AB.xml"), what = "", sep = "\n", quiet = TRUE, blank.lines.skip = FALSE)
jitter_value <- 0.001

clade_start <- grep("leaf.location24", xml_AB_rrw) 
jitter_start <- grep("</traitParameter>", xml_AB_rrw)
jitter_start_val <- jitter_start[jitter_start > clade_start][1]

xml_AB_rrw <- append(xml_AB_rrw, c(
    paste0('\t\t<jitter window=', '"', jitter_value, " ", jitter_value, '" duplicatesOnly="true">'),
    '\t\t\t<parameter idref="leaf.location24"/>',
    '\t\t</jitter>'
), after = jitter_start_val)

#writeLines(xml_AB_rrw, here(base_dir, "data", "BEAST", "EA-2021-AB", "RRW_run", "BEAST_RRW_EA-2021-AB.xml"))

## Run BEAST

genotype_vec <- c("EA-2020-A")
#genotype_vec <- c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB")

for(genotype in genotype_vec) {
    for(prefix in c("1_", "2_")) { # each XML-file is run twice
        local({
            wd <- getwd()
            on.exit(setwd(wd), add = TRUE)

            run_dir <- here(base_dir, "data", "BEAST", genotype, "RRW_run")
            xml <- paste0("BEAST_RRW_", genotype, ".xml")

            setwd(run_dir)

            system2(
                command = "beast",
                args = c("-overwrite", "-prefix", prefix ,xml), 
                stdout = file.path(run_dir, paste0(prefix, "beast_stdout", ".log")),
                stderr = file.path(run_dir, paste0(prefix, "beast_stderr", ".log")),
                wait = TRUE
            )
        })
        message("Finished BEAST RRW for genotype = ", genotype, ", prefix = ", prefix)
    }
}

## Combine log-files
burnin <- format(as.numeric(chain_length) / 100, scientific = FALSE, trim = TRUE)
#burnin <- "500000"

for (x in genotype_vec) {
    local({
        wd <- getwd()
        on.exit(setwd(wd), add = TRUE)

        run_dir <- here(base_dir, "data", "BEAST", x, "RRW_run")
        setwd(run_dir)

        ## 1) Combine LOG files
        ## logcombiner -burnin 5000000 1_All_clades.log 2_All_clades.log x_RRW.log 
        system2(
            command = "logcombiner",
            args = c(
                "-burnin", burnin, 
                "1_All_clades.log",
                "2_All_clades.log",
                paste0(x, "_RRW.log")
            ),
            stdout = file.path(run_dir, paste0(x, "_logcombiner_log_stdout.log")),
            stderr = file.path(run_dir, paste0(x, "_logcombiner_log_stderr.log")),
            wait = TRUE
        )
        message("Finished logcombiner (logs) for genotype ", x)

        ## 2) Combine TREE files
        ## logcombiner -trees -burnin 5000000 1_Clade_clade_index.trees 2_Clade_clade_index.trees Clade_clade_index.trees 
        clade_vec <- cluster_id_list[[x]]
        for(clade_index in clade_vec) {
            file_1 <- paste0("1_Clade_", clade_index, ".trees")
            file_2 <- paste0("2_Clade_", clade_index, ".trees")
            file_out <- paste0("Clade_", clade_index, ".trees")

            if (file.exists(file.path(run_dir, file_1)) & file.exists(file.path(run_dir, file_2))) {
                system2(
                    command = "logcombiner",
                    args = c(
                        "-trees",
                        "-burnin", burnin, 
                        paste0("1_Clade_", clade_index, ".trees"),
                        paste0("2_Clade_", clade_index, ".trees"),
                        paste0("Clade_", clade_index, ".trees")
                    ),
                    stdout = file.path(run_dir, paste0(x, "_logcombiner_trees_stdout.log")),
                    stderr = file.path(run_dir, paste0(x, "_logcombiner_trees_stderr.log")),
                    wait = TRUE
                )
            }
        }
        message("Finished logcombiner (trees) for genotype ", x)

        ## 3) Annotate trees using HIPSTR
        ## treeannotator -burnin 0 Clade_clade_index.trees Clade_clade_index.tree
        for(clade_index in clade_vec) {
            tree_file <- paste0("Clade_", clade_index, ".trees")
            out_file <- paste0("Clade_", clade_index, ".tree")

            if (file.exists(file.path(run_dir, tree_file))) {
                system2(
                    command = "treeannotator",
                    args = c(
                        "-burnin", "0",
                        "-hpd2D", "0.8",
                        paste0("Clade_", clade_index, ".trees"),
                        paste0("Clade_", clade_index, ".tree")
                    ),
                    stdout = file.path(run_dir, paste0(x, "_treeannotator_stdout.log")),
                    stderr = file.path(run_dir, paste0(x, "_treeannotator_stderr.log")),
                    wait = TRUE
                )
            }
        }
        message("Finished treeannotator for genotype ", x)
    })
}

## MCC instead of HIPSTR
### Note: MCC and HIPSTR trees are the same

for (x in genotype_vec) {
    local({
        wd <- getwd()
        on.exit(setwd(wd), add = TRUE)

        run_dir <- here(base_dir, "data", "BEAST", x, "RRW_run")
        setwd(run_dir)

        ## 3) Annotate trees using HIPSTR
        ## treeannotator -burnin 0 Clade_clade_index.trees Clade_clade_index.tree
        clade_vec <- cluster_id_list[[x]]

        for(clade_index in clade_vec) {
            tree_file <- paste0("Clade_", clade_index, ".trees")
            out_file <- paste0("Clade_", clade_index, ".tree")

            if (file.exists(file.path(run_dir, tree_file))) {
                system2(
                    command = "treeannotator",
                    args = c(
                        "-burnin", "0",
                        "-type", "mcc",
                        paste0("Clade_", clade_index, ".trees"),
                        paste0("Clade_", clade_index, "_mcc.tree")
                    ),
                    stdout = file.path(run_dir, paste0(x, "_treeannotator_mcc_stdout.log")),
                    stderr = file.path(run_dir, paste0(x, "_treeannotator_mcc_stderr.log")),
                    wait = TRUE
                )
            }
        }
        message("Finished treeannotator for genotype ", x)
    })
}

# Extracting spatio-temporal information embedded in HIPTSTR

MCC_tree_extractions = function(mcc_tre, mostRecentSamplingDatum) {
    mcc_tab = matrix(nrow=dim(mcc_tre$edge)[1], ncol=7)
    colnames(mcc_tab) = c("node1","node2","length","startLon","startLat","endLon","endLat")
    mcc_tab[,c("node1","node2")] = mcc_tre$edge
    mcc_tab[,c("length")] = mcc_tre$edge.length
    for (i in 1:length(mcc_tre$annotations)) {
        annotations = mcc_tre$annotations[[i]]
        mcc_tab[i,c("endLon","endLat")] = cbind(annotations$location.rate2, annotations$location.rate1)
    }
    for (i in 1:length(mcc_tre$annotations)) {
        index = which(mcc_tab[,"node2"] == mcc_tab[i,"node1"])
        if (length(index) > 0) {
                mcc_tab[i,c("startLon","startLat")] = mcc_tab[index,c("endLon","endLat")]
        } else {
                annotations = mcc_tre$root.annotation
                mcc_tab[i,c("startLon","startLat")] = cbind(annotations$location.rate2, annotations$location.rate1)
        }
    }
    l = length(mcc_tab[,1]); ll = matrix(1:l,nrow=l,ncol=l); ll[] = 0
    for (j in 1:l) {
        subMat = mcc_tab[j,2]
        subMat = subset(mcc_tab,mcc_tab[,2]==subMat)
        ll[j,1] = subMat[,3]
        subMat = subMat[1,1]
        subMat1 = subset(mcc_tab,mcc_tab[,2]==subMat)
        for (k in 1:l) {
            if (nrow(subMat1) > 0) {
                ll[j,k+1] = subMat1[,3]
                subMat2 = subMat1[1,1]
                subMat1 = subset(mcc_tab,mcc_tab[,2]==subMat2)
            }
        }
    }
    endNodeL = rowSums(ll)
    mcc_tab = cbind(mcc_tab, endNodeL)
    startNodeL = matrix(1:l,nrow=l,ncol=1)
    startNodeL[] = 0
    for (j in 1:l) {
        r = mcc_tab[j,1]
        s = subset(mcc_tab,mcc_tab[,2]==r)
        for (k in 1:l) {
            if (nrow(s) > 0) {
                startNodeL[j,1] = s[,8]
            }
        }	
    }
    colnames(startNodeL) = "startNodeL"
    mcc_tab = cbind(mcc_tab,startNodeL)
    maxEndLIndice = which.max(mcc_tab[,"endNodeL"])
    maxEndL = mcc_tab[maxEndLIndice,"endNodeL"]
    endYear = matrix(mcc_tab[,"endNodeL"]-maxEndL)
    endYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
    startYear = matrix(mcc_tab[,"startNodeL"]-maxEndL)
    startYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
    colnames(startYear) = "startYear"
    colnames(endYear) = "endYear"
    mcc_tab = cbind(mcc_tab,startYear,endYear)
    mcc_tab = mcc_tab[order(mcc_tab[,"startYear"],decreasing=F),]
    mcc_tab1 = mcc_tab[1,]
    mcc_tab2 = mcc_tab[2:dim(mcc_tab)[1],]
    mcc_tab2 = mcc_tab2[order(mcc_tab2[,"endYear"],decreasing=F),]
    mcc_tab = rbind(mcc_tab1, mcc_tab2)
    
    tipLabels = matrix(nrow = dim(mcc_tab)[1], ncol = 1)
    colnames(tipLabels) = "tipLabel"
    for (i in 1:length(mcc_tre$tip.label)) {
        tipLabels[which(mcc_tab[, "node2"] == i), "tipLabel"] = gsub("'", "", mcc_tre$tip.label[i])
    }
    mcc_tab = cbind(mcc_tab, tipLabels)

    return(mcc_tab)
}

genotype_vec  <- c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB")

for(genotype in genotype_vec) {
    run_dir <- here(base_dir, "data", "BEAST", genotype, "RRW_run")
    tree_files <- sub("\\.tree$", "", list.files(run_dir, pattern = "^Clade_[0-9]+\\.tree$", full.names = FALSE))

    all_clades <- NULL

    for(tree_file in tree_files) {
        clade_index <- as.numeric(sub("Clade_", "", tree_file))
        clade_data <- clusters_rrw[[genotype]][[which(sapply(clusters_rrw[[genotype]], function(x) x$cluster_id[1]) == clade_index)]]
        
        most_recent_sampling_datum <- max(as.numeric(clade_data$sampling_date))

        hipstr_tree <- seraphim::readAnnotatedNexus(file.path(run_dir, paste0(tree_file, ".tree")))
        rrw_output_tab <- MCC_tree_extractions(hipstr_tree, most_recent_sampling_datum)

        num_cols <- c("node1","node2","length","startLon","startLat",
              "endLon","endLat","endNodeL","startNodeL",
              "startYear","endYear")

        rrw_output_tab <- as.data.frame(rrw_output_tab)
        rrw_output_tab[num_cols] <- lapply(rrw_output_tab[num_cols], as.numeric)

        rrw_output_tab$cluster_id <- clade_index

        write.csv(rrw_output_tab, file.path(run_dir, paste0(tree_file, ".csv")), row.names = FALSE, quote = FALSE)

        if (is.null(all_clades)) {
            all_clades <- rrw_output_tab
        } else {
            maxNodeID <- max(all_clades[, c("node1", "node2")])
            rrw_output_tab[, c("node1", "node2")] <- rrw_output_tab[, c("node1", "node2")] + maxNodeID
            all_clades <- rbind(all_clades, rrw_output_tab)
        }
    }
    write.csv(all_clades, file.path(run_dir, "All_clades.csv"), row.names = FALSE, quote = FALSE)

    message("Finished MCC RRW extraction for ", genotype)
}

genotype_vec

for(genotype in genotype_vec) {
    run_dir <- here(base_dir, "data", "BEAST", genotype, "RRW_run")
    tree_files <- sub("\\.tree$", "", list.files(run_dir, pattern = "^Clade_[0-9]+\\.tree$", full.names = FALSE))

    if (length(tree_files) == 0) {
        message("Finished extracting information for ", genotype)
        next
    } # C2 has no files

    all_clades_ext_dir <- file.path(run_dir, "All_clades_ext")
    dir.create(all_clades_ext_dir, showWarnings = FALSE)

    trees_list <- setNames(vector("list", length(tree_files)), tree_files)
    most_recent_dates <- setNames(numeric(length(tree_files)), tree_files)
    n_trees_per_clade <- setNames(integer(length(tree_files)), tree_files)

    for (tree_file in tree_files) {
        clade_index <- as.numeric(sub("Clade_", "", tree_file))
        clade_data <- clusters_rrw[[genotype]][[which(sapply(clusters_rrw[[genotype]], function(x) x$cluster_id[1]) == clade_index)]]

        most_recent_dates[tree_file] <- max(as.numeric(clade_data$sampling_date), na.rm = TRUE)

        trees_list[[tree_file]] <- seraphim::readAnnotatedNexus(file.path(run_dir, paste0(tree_file, ".trees")))
        n_trees_per_clade[tree_file] <- length(trees_list[[tree_file]])
    }

    nberOfExtractionFiles <- min(n_trees_per_clade, na.rm = TRUE)

    for (i in seq_len(nberOfExtractionFiles)) {
        all_clades_post <- NULL

        for (tree_file in tree_files) {
            tab <- postTreeExtractions(post_tre = trees_list[[tree_file]][[i]], most_recent_dates[tree_file])

            tab[, "node1"] <- as.numeric(tab[, "node1"])
            tab[, "node2"] <- as.numeric(tab[, "node2"])

            if (is.null(all_clades_post)) {
                all_clades_post <- tab
            } else {
                maxNodeID <- max(c(as.numeric(all_clades_post[, "node1"]), as.numeric(all_clades_post[, "node2"])), na.rm = TRUE)

                tab[, "node1"] <- as.numeric(tab[, "node1"]) + maxNodeID
                tab[, "node2"] <- as.numeric(tab[, "node2"]) + maxNodeID

                all_clades_post <- rbind(all_clades_post, tab)
            }
        }

        write.csv(all_clades_post, file.path(all_clades_ext_dir, paste0("TreeExtractions_", i, ".csv")), row.names = FALSE, quote = FALSE)
    }
    message("Finished extracting information for ", genotype)
}

##################################################################################################################################################

# Plots
# Plotting spatio-temporal information in HIPSTR tree
## Estimate polygons/HPD regions associated with each internal node (instead of for each time slice)

genotype_vec <- c("EA-2020-A", "EA-2021-AB", "EA-2022-BB")
genotype_vec_short <- c("A", "AB", "BB")

annotation_edit_func <- function(x) { # location.rate instead of location
    nms <- names(x)
    nms <- sub("^location\\.rate1", "location1", nms)
    nms <- sub("^location\\.rate2", "location2", nms)
    names(x) <- nms
    x
}

## Global min and max date
global_minYear <- Inf
global_maxYear <- -Inf

for (genotype in genotype_vec) {
    run_dir <- here(base_dir, "data", "BEAST", genotype, "RRW_run")
    mcc_file <- file.path(run_dir, "All_clades.csv")

    lines <- readLines(mcc_file, warn = FALSE)
    if (length(lines) == 0 || all(trimws(lines) == "")) next

    mcc_tab <- read.csv(mcc_file)
    global_minYear <- min(global_minYear, mcc_tab$startYear, na.rm = TRUE)
    global_maxYear <- max(global_maxYear, mcc_tab$endYear, na.rm = TRUE)
}

## Colour scale
colour_scale <- colorRampPalette(
    c(
        "#5E3C99",  # purple
        "#2C7BB6",  # blue
        "#66C2A5",  # green
        "#FEE08B",  # darker yellow
        "#E66101"   # darker orange
    )
)(500)
n_col <- length(colour_scale)

## Loop for plotting
genotype_vec <- c("EA-2020-A", "EA-2021-AB", "EA-2022-BB")

for(i in seq_along(genotype_vec)) {
    genotype <- genotype_vec[i]
    genotype_s <- genotype_vec_short[i]

    message("Starting RRW plotting for ", genotype, " .....")

    run_dir <- here(base_dir, "data", "BEAST", genotype, "RRW_run")
    borders <- shapefile(here(base_dir, "plots", "gadm41_BEL_shp", "gadm41_BEL_2.shp"))

    localTreesDirectory <- file.path(run_dir, "All_clades_ext")
    
    mcc_tab <- read.csv(file.path(run_dir, "All_clades.csv"))

    names(mcc_tab)[names(mcc_tab) == "tipLabel"] <- "tip_id"
    tip_meta <- bind_rows(clusters_rrw[[genotype]]) %>% distinct(tip_id, species_group, cluster_id)

    mcc_tab <- left_join(mcc_tab, tip_meta, by = "tip_id")
    mcc_tab$point_pch <- c(21, 24, 23, 22, 13)[match(mcc_tab$species_group, c("wild_birds", "poultry", "mammals", "captive_birds", "other"))]
'    lty_vec <- c("solid", "dashed", "dotted", "dotdash", "longdash")[match(mcc_tab$cluster_id, unique(na.omit(mcc_tab$cluster_id)))]
    mcc_tab$lty <- NA
    for(n in seq_along(lty_vec)) {
        lty <- lty_vec[n]
        if(is.na(lty)) {
            next_non_na <- which(!is.na(lty_vec[(n+1):length(lty_vec)]))[1]
            if(!is.na(next_non_na)) mcc_tab$lty[n] <- lty_vec[n + next_non_na]
        } else mcc_tab$lty[n] <- lty
    }'

    polygons <- list()
    tree_file_list <- list.files(run_dir, pattern = "^Clade_[0-9]+\\.tree$", full.names = FALSE)
    
    for(clade_tree in tree_file_list) { # Loop across clades, get polygons for each clade, append in one list
        mcc_tree <- seraphim::readAnnotatedNexus(file.path(run_dir, clade_tree))
        
        mcc_tree$annotations <- lapply(mcc_tree$annotations, annotation_edit_func)
        mcc_tree$root.annotation <- annotation_edit_func(mcc_tree$root.annotation)

        polygons_clade <- suppressWarnings(spreadGraphic1(mcc_tree, mcc_tab)) # HPD regions associated with each internal node

        polygons <- c(polygons, polygons_clade)
    }

    nberOfExtractionFiles <- length(list.files(localTreesDirectory, pattern = "^TreeExtractions_[0-9]+\\.csv$"))
    startDatum <- min(mcc_tab$startYear, na.rm = TRUE)
    precision <- 1 / (365 / 7)
    cexNode <- 0.8
    percentage <- 80
    prob <- percentage / 100
    croppingPolygons <- FALSE
    
    #polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision)) # HPC regions per time slice

    minYear <- global_minYear
    maxYear <- global_maxYear

    startYears_indices <- round(((mcc_tab$startYear - minYear) / (maxYear - minYear)) * (n_col - 1)) + 1
    endYears_indices <- round(((mcc_tab$endYear - minYear) / (maxYear - minYear)) * (n_col - 1)) + 1

    startYears_indices[startYears_indices < 1] <- 1
    startYears_indices[startYears_indices > n_col] <- n_col
    endYears_indices[endYears_indices < 1] <- 1
    endYears_indices[endYears_indices > n_col] <- n_col

    startYears_colours <- colour_scale[startYears_indices]
    endYears_colours <- colour_scale[endYears_indices]

    polygons_colours <- rep(NA, length(polygons))
    for (j in seq_along(polygons)) {
        date <- as.numeric(names(polygons[[j]]))
        polygon_index <- round(((date - minYear) / (maxYear - minYear)) * (n_col - 1)) + 1
        polygon_index <- max(1, min(n_col, polygon_index))
        polygons_colours[j] <- adjustcolor(colour_scale[polygon_index], alpha.f = 0.25)
    }

    out_file <- here(base_dir, "plots", paste0("RRW_map_", genotype, ".png"))

    png(out_file, width = 6, height = 6.3, units = "in", res = 300)
    par(mar = c(0, 0, 0, 0), oma = c(0.5, 3.5, 1, 0), mgp = c(0, 0.4, 0), lwd = 0.2, bty = "o", bg = "grey99")

    plot(borders, col = "white", border = NA, axes = FALSE, xlim = c(bbox(borders)[1,1] - 0.3, bbox(borders)[1,2] + 0.3)) # white over grey background
    plot(raster::aggregate(borders), add = TRUE, col = NA, border = "grey10", lwd = 1) # dark grey international borders
    plot(borders, add = TRUE, lwd = 0.3, border = "gray30") # province borders

    for (j in length(polygons):1) {
        pol <- polygons[[j]]
        if (croppingPolygons) pol <- crop(pol, borders)
        plot(pol, add = TRUE, col = polygons_colours[j], border = NA, axes = FALSE)
    }

    for (j in seq_len(nrow(mcc_tab))) {
        curvedarrow(
            cbind(mcc_tab[j, "startLon"], mcc_tab[j, "startLat"]),
            cbind(mcc_tab[j, "endLon"], mcc_tab[j, "endLat"]),
            arr.length = 0, arr.width = 0, lwd = 0.5, #lty = mcc_tab$lty[j],
            lcol = "gray10", arr.col = NA, arr.pos = FALSE,
            curve = 0.1, dr = NA, endhead = FALSE
        )
    }

    for (j in rev(seq_len(nrow(mcc_tab)))) {
        if (!mcc_tab[j, "node1"] %in% mcc_tab[, "node2"]) {
            points(mcc_tab[j, "startLon"], mcc_tab[j, "startLat"],
                pch = 21, bg = startYears_colours[j], col = "gray10",
                cex = cexNode, lwd = 0.2)
        }
        points(mcc_tab[j, "endLon"], mcc_tab[j, "endLat"],
            pch = mcc_tab$point_pch[j],
            bg = endYears_colours[j], col = "gray10",
            cex = cexNode, lwd = 0.2)
    }

    legend_rast <- raster::raster(matrix(c(minYear, maxYear), nrow = 1))

    plot(
        legend_rast, legend.only = TRUE, add = TRUE, col = colour_scale,
        legend.width = 0.5, legend.shrink = 0.3, smallplot = c(0.08, 0.48, 0.14, 0.15),
        legend.args = list(text = "Year", cex = 0.7, line = 0.3, col = "gray30"),
        horizontal = TRUE,
        axis.args = list(cex.axis = 0.6, lwd = 0, lwd.tick = 0.2, tck = -0.8, col.axis = "gray30", line = 0, mgp = c(0, -0.02, 0), at = pretty(c(minYear, maxYear), n = 4))
    )

    legend("topright", inset = c(0.05, 0.05), legend = c("Wild birds", "Poultry", "Mammals", "Captive birds"), pch = c(21, 24, 23, 22), pt.bg = "grey50")

    dev.off()
    assign(paste0("rrw_plot_", genotype_s), out_file, envir = .GlobalEnv)
}

## Paneling
make_panel <- function(img_path, label) {
    img <- rasterGrob(readPNG(img_path), interpolate = TRUE)
    title <- textGrob(label, x = 0.04, y = 0.95, just = c("left", "top"), gp = gpar(fontsize = 14))  
    wrap_elements(full = grobTree(img, title))
}

p_A  <- make_panel(rrw_plot_A,  "EA-2020-A")
p_AB <- make_panel(rrw_plot_AB, "EA-2021-AB")
p_BB <- make_panel(rrw_plot_BB, "EA-2022-BB")

rrw_combined <- p_A | p_AB | p_BB

#ggsave(here(base_dir, "plots", "RRW_node_combined.png"), rrw_combined, width = 14, height = 5.2, dpi = 600, bg = "white")

# Plotting clusters individually per genotype
genotype_vec <- c("EA-2020-A", "EA-2021-AB", "EA-2022-BB")

rrw_cluster_panels <- list()

for(i in seq_along(genotype_vec)) {
    genotype <- genotype_vec[i]
    genotype_s <- genotype_vec_short[i]

    message("Starting RRW plotting for ", genotype, " .....")

    run_dir <- here(base_dir, "data", "BEAST", genotype, "RRW_run")
    borders <- shapefile(here(base_dir, "plots", "gadm41_BEL_shp", "gadm41_BEL_2.shp"))
    localTreesDirectory <- file.path(run_dir, "All_clades_ext")
    
    mcc_tab <- read.csv(file.path(run_dir, "All_clades.csv"))
    names(mcc_tab)[names(mcc_tab) == "tipLabel"] <- "tip_id"

    cexNode <- 0.8
    croppingPolygons <- FALSE
    minYear <- global_minYear
    maxYear <- global_maxYear

    tree_file_list <- list.files(run_dir, pattern = "^Clade_[0-9]+\\.tree$", full.names = FALSE)
    cluster_ids <- sort(as.numeric(gsub("^Clade_|\\.tree$", "", tree_file_list)))

    polygons_list <- list()
    for(clade_tree in tree_file_list) {
        clade_id <- as.numeric(gsub("^Clade_|\\.tree$", "", clade_tree))
        mcc_tree <- seraphim::readAnnotatedNexus(file.path(run_dir, clade_tree))
        
        mcc_tree$annotations <- lapply(mcc_tree$annotations, annotation_edit_func)
        mcc_tree$root.annotation <- annotation_edit_func(mcc_tree$root.annotation)

        polygons_list[[as.character(clade_id)]] <- suppressWarnings(spreadGraphic1(mcc_tree, mcc_tab))
    }

    cluster_plot_list <- vector("list", length(cluster_ids))

    for (k in seq_along(cluster_ids)) {
        cluster_id <- cluster_ids[k]

        mcc_tab_cluster <- read.csv(file.path(run_dir, paste0("Clade_", cluster_id, ".csv")))
        names(mcc_tab_cluster)[names(mcc_tab_cluster) == "tipLabel"] <- "tip_id"

        tip_meta_cluster <- bind_rows(clusters_rrw[[genotype]]) %>%
            distinct(tip_id, species_group, cluster_id) %>%
            filter(cluster_id == !!cluster_id)

        mcc_tab_cluster <- left_join(mcc_tab_cluster, tip_meta_cluster, by = "tip_id")
        mcc_tab_cluster$point_pch <- as.integer(c(21, 24, 23, 22, 13)[match(
            mcc_tab_cluster$species_group,
            c("wild_birds", "poultry", "mammals", "captive_birds", "other")
        )])

        startYears_indices_cluster <- round(((mcc_tab_cluster$startYear - minYear) / (maxYear - minYear)) * (n_col - 1)) + 1
        endYears_indices_cluster <- round(((mcc_tab_cluster$endYear - minYear) / (maxYear - minYear)) * (n_col - 1)) + 1

        startYears_indices_cluster[startYears_indices_cluster < 1] <- 1
        startYears_indices_cluster[startYears_indices_cluster > n_col] <- n_col
        endYears_indices_cluster[endYears_indices_cluster < 1] <- 1
        endYears_indices_cluster[endYears_indices_cluster > n_col] <- n_col

        mcc_tab_cluster$startYears_colours <- colour_scale[startYears_indices_cluster]
        mcc_tab_cluster$endYears_colours <- colour_scale[endYears_indices_cluster]

        polygons_cluster <- polygons_list[[as.character(cluster_id)]]

        polygons_colours <- rep(NA, length(polygons_cluster))
        for (j in seq_along(polygons_cluster)) {
            date <- as.numeric(names(polygons_cluster[[j]]))
            polygon_index <- round(((date - minYear) / (maxYear - minYear)) * (n_col - 1)) + 1
            polygon_index <- max(1, min(n_col, polygon_index))
            polygons_colours[j] <- adjustcolor(colour_scale[polygon_index], alpha.f = 0.25)
        }

        png_file <- tempfile(fileext = ".png")

        png(png_file, width = 6, height = 6.3, units = "in", res = 300)
        par(mar = c(0, 0, 0, 0), oma = c(0.5, 3.5, 1, 0), mgp = c(0, 0.4, 0), lwd = 0.2, bty = "o", bg = "grey99")

        plot(borders, col = "white", border = NA, axes = FALSE,
             xlim = c(bbox(borders)[1,1] - 0.3, bbox(borders)[1,2] + 0.3))
        plot(raster::aggregate(borders), add = TRUE, col = NA, border = "grey10", lwd = 1)
        plot(borders, add = TRUE, lwd = 0.3, border = "gray30")

        for (j in seq_along(polygons_cluster)) {
            pol <- polygons_cluster[[j]]
            if (croppingPolygons) pol <- crop(pol, borders)
            plot(pol, add = TRUE, col = polygons_colours[j], border = NA, axes = FALSE)
        }

        for (j in seq_len(nrow(mcc_tab_cluster))) {
            if (any(is.na(mcc_tab_cluster[j, c("startLon","startLat","endLon","endLat")]))) next

            curvedarrow(
                cbind(mcc_tab_cluster[j, "startLon"], mcc_tab_cluster[j, "startLat"]),
                cbind(mcc_tab_cluster[j, "endLon"], mcc_tab_cluster[j, "endLat"]),
                arr.length = 0, arr.width = 0, lwd = 0.5,
                lcol = "gray10", arr.col = NA, arr.pos = FALSE,
                curve = 0.1, dr = NA, endhead = FALSE
            )
        }

        for (j in rev(seq_len(nrow(mcc_tab_cluster)))) {
            if (!mcc_tab_cluster[j, "node1"] %in% mcc_tab_cluster[, "node2"] &&
                !any(is.na(mcc_tab_cluster[j, c("startLon","startLat")]))) {
                points(mcc_tab_cluster[j, "startLon"], mcc_tab_cluster[j, "startLat"],
                       pch = 21, bg = mcc_tab_cluster$startYears_colours[j], col = "gray10",
                       cex = cexNode, lwd = 0.2)
            }

            if (!any(is.na(mcc_tab_cluster[j, c("endLon","endLat")]))) {
                points(mcc_tab_cluster[j, "endLon"], mcc_tab_cluster[j, "endLat"],
                       pch = mcc_tab_cluster$point_pch[j], bg = mcc_tab_cluster$endYears_colours[j],
                       col = "gray10", cex = cexNode, lwd = 0.2)
            }
        }

        legend_rast <- raster::raster(matrix(c(minYear, maxYear), nrow = 1))

        plot(
            legend_rast, legend.only = TRUE, add = TRUE, col = colour_scale,
            legend.width = 0.5, legend.shrink = 0.3, smallplot = c(0.08, 0.48, 0.14, 0.15),
            legend.args = list(text = "Year", cex = 0.7, line = 0.3, col = "gray30"),
            horizontal = TRUE,
            axis.args = list(cex.axis = 0.6, lwd = 0, lwd.tick = 0.2, tck = -0.8, col.axis = "gray30", line = 0, mgp = c(0, -0.02, 0), at = pretty(c(minYear, maxYear), n = 4))
        )

        legend("topright", inset = c(0.05, 0.05),
               legend = c("Wild birds", "Poultry", "Mammals", "Captive birds"),
               pch = c(21, 24, 23, 22), pt.bg = "grey50")

        dev.off()

        img <- rasterGrob(readPNG(png_file), interpolate = TRUE)
        cluster_plot_list[[k]] <- wrap_elements(full = img)
    }

    rrw_cluster_panels[[genotype]] <- wrap_plots(cluster_plot_list)

}

#ggsave(here(base_dir, "plots", "RRW_clusters_EA-2021-AB.png"),
#       rrw_cluster_panels[["EA-2021-AB"]], width = 12, height = 8, dpi = 600, bg = "white")

# Dispersal statistics based on the continuous phylogeographic reconstructions

for(genotype in genotype_vec) {
    message("Estimating dispersal statistics 1 for ", genotype, " .....")

    run_dir <- here(base_dir, "data")
    localTreesDirectory <- file.path(run_dir, "BEAST", genotype, "RRW_run", "All_clades_ext")
    nberOfExtractionFiles <- length(list.files(localTreesDirectory, pattern = "^TreeExtractions_[0-9]+\\.csv$"))

    if(nberOfExtractionFiles == 0) next

    timeSlices <- 100; onlyTipBranches <- FALSE; showingPlots <- FALSE; nberOfCores <- 5; slidingWindow  <- 1
    outputName <- paste0(run_dir, "/RRW_stat/", genotype, "/", genotype)

    seraphim::spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
}

rrw_stat1_df <- read.table(here(base_dir, "data", "RRW_stat", "RRW_stat.txt"), header = TRUE)

rrw_stat2_list <- list()

for(genotype in genotype_vec) {
    run_dir <- here(base_dir, "data", "RRW_stat", genotype)
    stats_file <- here(run_dir, paste0(genotype, "_estimated_dispersal_statistics.txt"))

    if (!file.exists(stats_file)) next

    stats <- read.table(stats_file, header = TRUE)

    wldv <- stats[,"weighted_branch_dispersal_velocity"] / 365
    wldv_md <- round(median(wldv), 2)
    wldv_qs <- round(quantile(wldv, c(0.025, 0.975)), 2)
    wldv_hpd <- round(HDInterval::hdi(wldv)[1:2], 2)

    wdc <- stats[,"weighted_diffusion_coefficient"]
    wdc_md <- round(median(wdc), 1)
    wdc_qs <- round(quantile(wdc, c(0.025, 0.975)), 1)
    wdc_hpd <- round(HDInterval::hdi(wdc)[1:2], 2)

    cat(genotype, ":\tWLDV = ", wldv_md, " km/day [", wldv_hpd[1], "-", wldv_hpd[2], "]\t\tWDC = ", wdc_md, " km2/year [", wdc_hpd[1], "-", wdc_hpd[2], "]\n", sep = "")

    rrw_stat2_list[[genotype]] <- data.frame(
        genotype = genotype,
        WLDV_median = wldv_md,
        WLDV_HPD_lower = wldv_hpd[1],
        WLDV_HPD_upper = wldv_hpd[2],
        WDC_median = wdc_md,
        WDC_HPD_lower = wdc_hpd[1],
        WDC_HPD_upper = wdc_hpd[2]
    )
}

'
EA-2020-A:      WLDV = 0.44 km/day [0.35-0.68]          WDC = 5186 km2/year [3241.64-10710.22]
EA-2021-AB:     WLDV = 0.71 km/day [0.62-0.91]          WDC = 5864.3 km2/year [4686.49-8089.78]
EA-2022-BB:     WLDV = 0.84 km/day [0.59-1.42]          WDC = 2682.6 km2/year [1272.02-7851.4]
'

rrw_stat2_df <- do.call(rbind, rrw_stat2_list)

# Introduction dynamics

## Histogram with clade size
head(cluster_size_table)

cluster_size_table_edit <- cluster_size_table
cluster_size_table_edit$genotype <- ifelse(cluster_size_table_edit$genotype == "EA-2020-C2", "EA-2020-C", cluster_size_table_edit$genotype)

histo_cluster_size <- ggplot(cluster_size_table_edit, aes(x = n_sequences, fill = genotype)) +
    geom_histogram(binwidth = 1, boundary = 0.5, colour = "black") +
    guides(fill = guide_legend(position = "inside")) +
    scale_x_continuous(breaks = seq(1, 20, by = 1)) +
    scale_y_continuous(breaks = seq(0, 150, by = 10)) +
    scale_fill_manual(values = met.brewer("Hokusai3", direction = 1)) +
    theme_classic() +
    theme(
        legend.position.inside = c(0.9, 0.85),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    )

#ggsave(here(base_dir, "plots", "histo_cluster_size.png"), width = 25, height = 12, units = "cm", dpi = 600)

## Histogram with introduction time of clades
colour_scale

clade_data_df <- bind_rows(
  lapply(names(clusters2), function(genotype) {
    bind_rows(
      lapply(clusters2[[genotype]], function(clade) {
        df <- bind_rows(clade) %>%
          filter(!is.na(sampling_date))
        
        if (nrow(df) == 0) return(NULL)
        
        df$genotype <- genotype
        df$clade_id <- unique(df$cluster_id)
        df
      }) 
    )
  })
)
rownames(clade_data_df) <- NULL

all_intro_date_df <- data.frame()

genotype_vec <- c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB")

node_loc_list <- list()

for (genotype in genotype_vec) {
    run_dir <- here(base_dir, "data", "BEAST", genotype, "DTA_run")
    
    beast_tree <- treeio::read.beast(file.path(run_dir, paste0(genotype, "_DTA.tree")))
    tree <- beast_tree@phylo
    tree_annotation <- as.data.frame(beast_tree@data)
    
    # get edges
    tree_edge <- tree$edge
    tree_parent_node <- tree_edge[,1]
    tree_child_node <- tree_edge[,2]
    
    # get location states
    ## node_loc_list for nodes for each genotype (see later)
    node_loc_list[[genotype]] <- setNames(as.character(tree_annotation$location.rate), as.character(tree_annotation$node))    
    start_location <- node_loc_list[[genotype]][as.character(tree_parent_node)]
    end_location <- node_loc_list[[genotype]][as.character(tree_child_node)]
    
    # define introduction edges
    intro_edges <- which(end_location == "Belgium" & start_location == "other")
    if (length(intro_edges) == 0) next
    
    # compute depths
    depth <- ape::node.depth.edgelength(tree)
    
    intro_parent <- tree$edge[intro_edges, 1]
    intro_child <- tree$edge[intro_edges, 2]
    intro_bl <- tree$edge.length[intro_edges]
    intro_mid_depth <- depth[intro_parent] + intro_bl / 2
    
    # tip dates
    tip_meta <- clade_data_df[clade_data_df$genotype == genotype, c("tip_id", "sampling_date")]
    tip_meta <- tip_meta[!duplicated(tip_meta$tip_id), ]
    
    tip_dates <- tip_meta$sampling_date
    names(tip_dates) <- tip_meta$tip_id
    tip_dates <- tip_dates[tree$tip.label]
    
    # root date and intro dates
    root_date <- median(tip_dates - depth[1:ape::Ntip(tree)], na.rm = TRUE)
    intro_date <- root_date + intro_mid_depth
    
    intro_date_df <- data.frame(
        genotype = genotype,
        edge_row = intro_edges,
        parent_node = intro_parent,
        child_node = intro_child,
        intro_date = intro_date
    )
    
    all_intro_date_df <- rbind(all_intro_date_df, intro_date_df)
}
rownames(all_intro_date_df) <- NULL

all_intro_date_df$date <- date_decimal(all_intro_date_df$intro_date)
all_intro_date_df$genotype <- ifelse(all_intro_date_df$genotype == "EA-2020-C2", "EA-2020-C", all_intro_date_df$genotype)

histo_intro_date <- ggplot(all_intro_date_df, aes(x = intro_date, fill = genotype)) +
    geom_histogram(
        binwidth = 1/12,   # one month in decimal years
        center = 2020 + 4/12 + 1/24,
        colour = "black"
    ) +
    scale_x_continuous(
        breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 2/12),
        minor_breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 1/12),
        labels = function(x) format(lubridate::date_decimal(x), "%m-%Y")
    ) +
    scale_y_continuous(breaks = seq(0, 21, by = 2)) +
    scale_fill_manual(values = met.brewer("Hokusai3", direction = 1)) +
    guides(
        x = guide_axis(minor.ticks = TRUE),
        fill = guide_legend(position = "inside")
    ) +
    theme_classic() +
    theme(
        legend.position.inside = c(0.9, 0.85),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

#ggsave(here(base_dir, "plots", "histo_intro_date.png"), width = 25, height = 12, units = "cm", dpi = 600)

ggplot(all_intro_date_df, aes(x = intro_date)) +
    geom_histogram(
        aes(fill = after_stat(x)),
        binwidth = 1/12,
        center = 2020 + 4/12 + 1/24,
        colour = "black"
    ) +
    guides(fill = "none") +
    scale_x_continuous(
        breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 1/12),
        labels = function(x) format(lubridate::date_decimal(x), "%m-%Y")
    ) +
    scale_y_continuous(breaks = seq(0, 150, by = 1)) +
    scale_fill_gradientn(colours = colour_scale) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
    
'ggsave(here(base_dir, "plots", "histo_intro_date_colour.png"), width = 25, height = 12, units = "cm", dpi = 600)

ggsave(here(base_dir, "plots", "histo_panel_2.png"), 
    histo_cluster_size + histo_intro_date +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold")),
    width = 12, height = 4, dpi = 600, units = "in"
)'

# Cluster visualisation: Tip nodes and centroid 
## Put all tip nodes on map and connect to centroid with curved dashed line to show spread around centroid for one cluster
## As before: All tips coloured according to collection date, different symbol for poultry and wild birds 
## Centroid: x = mean(x_all); y = mean(y_all)

genotype_vec <- c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB")

shapefile <- read_sf(here(base_dir, "data", "metadata", "gadm41_BEL_2.shp"))

plot_colours <- met.brewer("Lakota", direction = -1)

tab_nodes_list <- lapply(genotype_vec, function(genotype) {
    bind_rows(
    lapply(seq_along(clusters2[[genotype]]), function(i) {
        clusters2[[genotype]][[i]] %>% mutate(genotype = genotype, cluster_id = i)
    })
    ) %>%
    filter(!is.na(latitude)) %>%
    dplyr::select(genotype, cluster_id, tip_id, sampling_date, latitude, longitude, species_group) %>%
    mutate(species_group = factor(species_group, levels = c("poultry", "wild_birds", "mammals", "captive_birds"))) %>%    add_count(cluster_id) %>%
    mutate(cluster_id2 = case_when(
        n == 1 ~ -1,
        n == 2 ~ 0,
        n > 2 ~ cluster_id
    )) %>%
    group_by(cluster_id) %>%
    mutate(
        mean_lon = if_else(n[1] >= 2, mean(longitude), NA_real_),
        mean_lat = if_else(n[1] >= 2, mean(latitude), NA_real_)
    ) %>%
    ungroup() %>%
    dplyr::select(-n) %>%
    group_by(cluster_id, latitude, longitude) %>%
    filter(!(cluster_id2 == 0 & row_number() > 1)) %>%
    ungroup() %>%
    group_by(cluster_id) %>%
    filter(!(cluster_id2 == 0 & n() == 1)) %>%
    mutate(cluster_id2 = if_else(cluster_id2 == 0 & n() == 1, -1L, cluster_id2)) %>%
    ungroup()
}) %>% setNames(genotype_vec)

tab_centroids <- lapply(genotype_vec, function(genotype) {
    tab_nodes_list[[genotype]] %>%
    distinct(cluster_id, cluster_id2, mean_lon, mean_lat)
}) %>% setNames(genotype_vec)

print(tab_nodes_list[[3]], n = 50)

plot_list <- lapply(genotype_vec, function(genotype) {  
    ggplot() +
        geom_sf(data = shapefile, fill = "white", colour = "grey50", linewidth = 0.2) +
        geom_sf(data = st_union(shapefile), fill = NA, colour = "grey20", linewidth = 0.3) +
        {if (nrow(tab_nodes_list[[genotype]] %>% filter(!is.na(mean_lon))) > 0)
        # lines
            geom_curve(data = tab_nodes_list[[genotype]] %>% filter(!is.na(mean_lon) & cluster_id2 >= 0),
                aes(x = longitude, y = latitude, xend = mean_lon, yend = mean_lat,
                colour = factor(cluster_id2)), curvature = 0.0, linetype = "dotted", alpha = 0.5,
                show.legend = FALSE)
        } +
        # data points
        geom_point(data = tab_nodes_list[[genotype]], 
            aes(x = longitude, y = latitude, fill = factor(cluster_id2), colour = factor(cluster_id2), shape = species_group), 
            size = 2.5, show.legend = TRUE) + 
        # centroids
        geom_point(data = tab_centroids[[genotype]], 
            aes(x = mean_lon, y = mean_lat, colour = factor(cluster_id2)), 
             size = if_else(tab_centroids[[genotype]]$cluster_id2 == 0, 1, 2), 
             stroke = if_else(tab_centroids[[genotype]]$cluster_id2 == 0, 1, 1.5), shape = 8, show.legend = FALSE) +
        scale_colour_manual(
            values = setNames(
                c("grey65", "grey65", plot_colours),
                c(-1, 0, sort(unique(tab_nodes_list[[genotype]]$cluster_id2[tab_nodes_list[[genotype]]$cluster_id2 > 0])))
            ),
            guide = "none",
            aesthetics = c("colour", "fill")
        ) +  
        scale_shape_manual(
            values = c("poultry" = 21, "wild_birds" = 24, "mammals" = 23, "captive_birds" = 22),
            breaks = c("poultry", "wild_birds", "mammals", "captive_birds"),
            labels = c("Poultry", "Wild birds", "Mammals", "Captive birds"),
            name = "Species group",
            drop = FALSE
        ) +
        guides(
            fill = "none",
            colour = "none",
            shape = guide_legend(override.aes = list(colour = "black", size = 4))
        )+
#        facet_wrap(~ cluster_id) +
        labs(title = ifelse(genotype == "EA-2020-C2", "EA-2020-C", genotype)) +
        theme_bw() +
        theme(
            panel.background = element_rect(fill = "grey95", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA),
            legend.key = element_rect(fill = NA, colour = NA),
            strip.background = element_rect(fill = "white", colour = NA),
            legend.position = "top",
            plot.title = element_text(face = "bold")
    )
})

i <- 3
plot_list[i]
genotype <- genotype_vec[i]
#ggsave(here(base_dir, "plots", paste0(genotype, "_centroid.png")), plot_list[i], width = 12, height = 8, dpi = 600, bg = "white")

plot_list[2] + labs(title = "EA-2020-C")

'ggsave(here(base_dir, "plots", "centroid_panel.png"), 
    wrap_plots(plot_list, ncol = 2) +
    plot_annotation(tag_levels = list(names(plot_list))) &
    theme(plot.tag = element_text(face = "bold")),
    width = 14, height = 12, dpi = 600, units = "in"
)'

'ggsave(here(base_dir, "plots", "histo_centroid_panel.png"),
  ((histo_cluster_size | histo_intro_date) / wrap_plots(plot_list, ncol = 2)) +
    plot_layout(height = c(1, 3)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold")),
  width = 12, height = 15, dpi = 300, units = "in"
)'

# Percentage of species per genotype
species_pct <- dplyr::bind_rows(lapply(names(clusters2), function(g) {
  if (length(clusters2[[g]]) == 0) return(NULL)
  dplyr::bind_rows(clusters2[[g]]) %>%
    filter(!is.na(species_group)) %>%
    dplyr::count(species_group) %>%
    dplyr::mutate(
      genotype = g,
      percentage = 100 * n / sum(n)
    )
}))

# Histogram: Epi data
epi_data

species_colours_1 <- c(
    "wild_birds" = met.brewer("Hokusai1")[5],
    "poultry" = met.brewer("Hokusai1")[4],
    "captive_birds" = met.brewer("Hokusai1")[6],
    "no_bird" = met.brewer("Hokusai1")[7]
)

histo_epi <- ggplot(epi_data, aes(x = decimal_date, fill = species_group)) +
    geom_histogram(
        binwidth = 1/12,
        center = 2020 + 4/12 + 1/24,
        colour = "black"
    ) +
    scale_x_continuous(
        breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 2/12),
        minor_breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 1/12),
        labels = function(x) format(lubridate::date_decimal(x), "%m-%Y")
    ) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    scale_fill_manual(
        values = species_colours_1,
        breaks = c("wild_birds", "poultry", "captive_birds", "no_bird"),
        labels = c("Wild birds", "Poultry", "Captive birds", "Not bird"),
        drop = FALSE
    ) +
    guides(
        x = guide_axis(minor.ticks = TRUE),
        fill = guide_legend(position = "inside")
    ) +
    theme_classic() +
    theme(
        legend.position.inside = c(0.9, 0.85),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

#ggsave(here("plots", "histo_epi.png"), width = 25, height = 12, units = "cm", dpi = 600)

species_colours_2 <- c(
    "wild_birds" = met.brewer("Hokusai1")[5],
    "poultry" = met.brewer("Hokusai1")[4],
    "captive_birds" = met.brewer("Hokusai1")[6],
    "mammals" = met.brewer("Hokusai1")[7]
)

histo_seq <- bind_rows(lapply(genotype_vec, function(g) {
  do.call(bind_rows, clusters2[[g]]) %>%
    filter(!is.na(species_group)) %>%
    mutate(species_group = factor(species_group, levels = c("wild_birds", "poultry", "captive_birds", "mammals")))
})) %>%
  ggplot(aes(x = sampling_date, fill = species_group)) +
  geom_histogram(
    binwidth = 1/12,
    center = 2020 + 4/12 + 1/24,
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 2/12),
    minor_breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 1/12),
    labels = function(x) format(lubridate::date_decimal(x), "%m-%Y")
  ) +
  scale_y_continuous(breaks = seq(0, 21, by = 2)) +
  scale_fill_manual(
        values = species_colours_2,
        breaks = c("wild_birds", "poultry", "captive_birds", "mammals"),
        labels = c("Wild birds", "Poultry", "Captive birds", "Mammals"),
        drop = FALSE
  ) +
  guides(
    x = guide_axis(minor.ticks = TRUE),
    fill = guide_legend(position = "inside")
  ) +
  theme_classic() +
  theme(
    legend.position.inside = c(0.9, 0.85),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

histo_seq_genotype <- bind_rows(lapply(genotype_vec, function(g) {
  do.call(bind_rows, clusters2[[g]]) %>%
    filter(!is.na(latitude)) %>%
    mutate(
      genotype = factor(
        g,
        levels = c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB")
      )
    )
})) %>%
  ggplot(aes(x = sampling_date, fill = genotype)) +
  geom_histogram(
    binwidth = 1/12,
    center = 2020 + 4/12 + 1/24,
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 2/12),
    minor_breaks = seq(2020 + 4/12 + 1/24, 2023 + 10/12 + 1/24, by = 1/12),
    labels = function(x) format(lubridate::date_decimal(x), "%m-%Y")
  ) +
  scale_y_continuous(breaks = seq(0, 21, by = 2)) +
  scale_fill_manual(
    values = MetBrewer::met.brewer("Hokusai3"),
    breaks = c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB"),
    labels = c("EA-2020-A", "EA-2020-C", "EA-2021-AB", "EA-2022-BB"),
    drop = FALSE
  ) +
  guides(
    x = guide_axis(minor.ticks = TRUE),
    fill = guide_legend(position = "inside")
  ) +
  theme_classic() +
  theme(
    legend.position.inside = c(0.9, 0.85),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

#ggsave(here("plots", "histo_seq_genotype.png"), width = 25, height = 12, units = "cm", dpi = 600)

# Percentage of genotype over time

bind_rows(lapply(genotype_vec, function(g) {
  do.call(bind_rows, clusters2[[g]]) %>%
    mutate(genotype = g)
})) %>%
  mutate(
    month = floor_date(date_decimal(sampling_date), "month"),
    genotype = factor(
      genotype,
      levels = c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB")
    )
  ) %>%
  dplyr::count(month, genotype) %>%
  ggplot(aes(month, n, fill = genotype)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = met.brewer("Hokusai3"),
    breaks = c("EA-2020-A", "EA-2020-C2", "EA-2021-AB", "EA-2022-BB"),
    labels = c("EA-2020-A", "EA-2020-C", "EA-2021-AB", "EA-2022-BB"),
    drop = FALSE
  ) +
  scale_x_date(
  breaks = seq(as.Date("2020-05-01"), as.Date("2023-12-01"), by = "2 months"),
  minor_breaks = seq(as.Date("2020-05-01"), as.Date("2023-12-01"), by = "1 month"),
  labels = function(x) format(x, "%m-%Y")
    ) +
    guides(
    x = guide_axis(minor.ticks = TRUE)
    ) +
    theme_classic() +
    theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
    )
)
#ggsave(here("plots", "percent_seq_genotype.png"), width = 25, height = 12, units = "cm", dpi = 600)

# Paneling
x_limits <- c(2020 + 5/12, 2024)

histo_epi <- histo_epi + coord_cartesian(xlim = x_limits)
histo_seq <- histo_seq + coord_cartesian(xlim = x_limits)

histo_panel <- histo_epi / histo_seq / histo_intro_date / histo_cluster_size
histo_panel_2 <- (histo_epi / histo_seq | histo_intro_date / histo_cluster_size )+ plot_annotation(tag_levels = "A")
#ggsave(here("plots", "histo_panel_2.png"), width = 40, height = 20, units = "cm", dpi = 600)

map_panel <- wrap_plots(plot_list, ncol = 2)

'ggsave(here(base_dir, "plots", "histo_centroid_panel_2.pdf"),
  (histo_panel | map_panel) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold")),
  width = 58, height = 35, dpi = 300, units = "cm"
)
'

# Calculate correlation between outbreak data and sequence data 

## Outbreak data: count reported outbreaks per month
epi_monthly <- epi_data %>%
  mutate(
    date = dmy(final_date),
    month = floor_date(date, unit = "month")
  ) %>%
  dplyr::count(month, name = "n_outbreaks")

## Sequence data: combine all genotypes and clusters
seq_data <- clusters2 %>%
  imap_dfr(~ bind_rows(.x) %>%
             mutate(genotype = .y))

## Keep Belgian sequences with coordinates and count per month
seq_monthly <- seq_data %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  mutate(
    date = as.Date(date_decimal(sampling_date)),
    month = floor_date(date, unit = "month")
  ) %>%
  dplyr::count(month, name = "n_sequences")

## Combine outbreak and sequence counts
correlation_data <- full_join(epi_monthly, seq_monthly, by = "month") %>%
  complete(
    month = seq(
      min(month, na.rm = TRUE),
      max(month, na.rm = TRUE),
      by = "month"
    ),
    fill = list(
      n_outbreaks = 0,
      n_sequences = 0
    )
  ) %>%
  mutate(
    n_outbreaks = replace_na(n_outbreaks, 0),
    n_sequences = replace_na(n_sequences, 0)
  ) %>%
  arrange(month)

correlation_model <- lm(
  n_sequences ~ n_outbreaks,
  data = correlation_data
)

summary(correlation_model)$r.squared

ggplot(correlation_data, aes(x = n_outbreaks, y = n_sequences)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Number of reported outbreaks per month",
    y = "Number of sequences per month"
  ) +
  theme_classic()

#ggsave(here("plots", "seq_epi_correlation.png"), width = 20, height = 14, units = "cm", dpi = 600)

# Pie charts for species groups

species_BE_df <- bind_rows(clusters2) %>%
    dplyr::count(species_group, species_host) %>% 
    as.data.frame()

species_BE_genotype_df <- bind_rows(
    lapply(names(clusters2), function(genotype) {
        bind_rows(clusters2[[genotype]]) %>%
            dplyr::count(species_group, species_host) %>%
            mutate(genotype = genotype)
    })
) %>% dplyr::select(genotype, species_group, species_host, n)

species_BE_table <- as.data.frame.matrix(table(species_BE_genotype_df$genotype, species_BE_genotype_df$species_group)) %>% 
    mutate(genotype = row.names(.)) %>% dplyr::relocate(genotype)

species_BE_table %>%
    pivot_longer(-genotype, names_to = "species_group", values_to = "n") %>%
    ggplot(aes(x = 1, y = n, fill = species_group)) +
    geom_col() +
    coord_polar(theta = "y") +
    facet_wrap(~genotype) +
    theme_void()

species_BE_genotype_df %>%
    drop_na(species_group) %>%
    dplyr::count(genotype, species_group, wt = n) %>%
    group_by(genotype) %>%
    mutate(pct = n / sum(n)) %>%
    ggplot(aes(1, pct, fill = species_group)) +
    geom_col(colour = "white") +
    geom_text(aes(label = scales::percent(pct, 0.1)),
        position = position_stack(vjust = 0.5), size = 3) +
    coord_polar("y") +
    facet_wrap(~ genotype) +
    scale_fill_manual(values = species_colours_2) +
    theme_void()
