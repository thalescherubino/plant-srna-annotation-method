"""Species profiles for GEO library curation and selection (curate_geo_candidates.py, select_libraries.py).

Each profile defines, as explicit and reviewable regular expressions:
  target      organism name that must appear in the GEO organism field (anything else is non_target);
  genotype()  the "subspecies" column: cultivar/genotype group of the sample, or non_target;
  rank        selection preference among those groups (lower first; reference genotype first);
  organs      ordered (category, pattern) list; the first match wins;
  deleterious obviously deleterious lines (sterile etc.), always excluded;
  extra_mutant / extra_other_line  species-specific mutant or special-line terms added to the shared rules.
Rice (2026-11-06) is unchanged from the original scripts. Tomato added 2026-11-09.
"""
import re


def _rice_genotype(r, text):
    org = r['organism'].lower()
    if 'oryza sativa' not in org:
        return 'non_target'
    if re.search(r'nipp?on?bare|nipponbare|nippnbare', text):
        return 'japonica_nipponbare'
    if 'indica' in org or re.search(r'indica|93-?11|9311|ir24|ir64|irbb|pusa basmati|pb-?1\b|zhenshan|minghui|teqing|kasalath|co39|n22|aus\b', text):
        return 'indica_or_aus'
    if 'japonica' in org or re.search(r'japonica|japanica|japponica|kitaake|zhonghua|dongjin|hwayoung|nongken|taichung|'
                                      r'koshihikari|wuyunjing|ilmi|illmi|jowoon|gopum|azucena|\baz_|9522|zhonghua ?11|'
                                      r'heukgyeong|nipponbare', text):
        return 'japonica_other'
    return 'unknown_sativa'


RICE_ORGANS = [
    ('mixed_tissues', r',.*\b(and|,)\b.*|inflorescences, leaves|leaves, roots|mixed|pooled tissues'),
    ('anther_pollen_meiocyte', r'anther|pollen|meiocyte|microspore|tetrad|stamen'),
    ('root', r'\broots?\b|radicle|root tip'),
    ('panicle_inflorescence', r'panicle|inflorescence|spikelet|floret|flower|young ear|rachis|glume|pistil|ovar'),
    ('seed_embryo_endosperm', r'embryo|endosperm|\bseeds?\b|grain|caryopsis|\bdap\b|aleurone'),
    ('callus_culture', r'callus|cell culture|suspension|protoplast'),
    ('leaf', r'flag ?leaf|leaf|leaves|blade'),
    ('sheath_stem_node', r'sheath|stem|culm|node|internode|tiller|crown|shoot base|shoot apex|\bsam\b|meristem'),
    ('seedling_shoot', r'seedlings?|shoot|aerial|whole plant|plantlet|coleoptile'),
]

# Tomato: other Solanum species (and potato) are non_target unless the sample is an introgression line or an
# interspecific hybrid in an S. lycopersicum background, which is a special line (interspecific_line).
TOMATO_WILD = re.compile(r'pennellii|pimpinellifolium|habrochaites|hirsutum|chilense|peruvianum|lycopersicoides|'
                         r'chmielewskii|neorickii|arcanum|corneliomulleri|galapagense|cheesmaniae|tuberosum|'
                         r'nicotiana|benthamiana', re.I)
TOMATO_IL = re.compile(r'introgression|\bil ?\d|\bils?\b|\bbil\b|hybrid|\bf1\b|x s\. ?|× ?s\.', re.I)


def _tomato_genotype(r, text):
    org = r['organism'].lower()
    if 'solanum lycopersicum' not in org:
        return 'non_target'
    if TOMATO_WILD.search(text):
        return 'interspecific_line' if TOMATO_IL.search(text) else 'non_target'
    if re.search(r'heinz|\bh ?1706\b', text):
        return 'heinz_1706'
    if re.search(r'cerasiforme|cherry', text) or 'cerasiforme' in org:
        return 'cerasiforme'
    if re.search(r'\bm ?82\b|micro-? ?tom|microtom|ailsa ?craig|\bac\b|money ?maker|castlemart|rio ?grande|'
                 r'moneyberg|alisa ?craig|zhongshu|\bvf36\b|san ?marzano|marmande|\bla\d{3,4}\b|\btr\d|cultivar', text):
        return 'cultivar_other'
    return 'unknown_lycopersicum'


TOMATO_ORGANS = [
    ('mixed_tissues', r',.*\b(and|,)\b.*|mixed|pooled tissues'),
    ('anther_pollen_meiocyte', r'anther|pollen|meiocyte|microspore|tetrad|stamen'),
    ('root', r'\broots?\b|radicle|root tip|hairy root'),
    ('fruit', r'fruit|pericarp|mature green|immature green|breaker|red ripe|ripening|turning stage|pink stage|'
              r'locular|placenta|columella|pulp|edible portion|\bdpa\b|\bmg\b|\bbr\b|\brr\b'),
    ('ovary_ovule', r'ovar|ovule|pistil|carpel|gynoecium'),
    ('flower_inflorescence', r'flower|inflorescence|floral|petal|sepal|truss|\bbuds?\b'),
    ('seed_embryo_endosperm', r'embryo|endosperm|\bseeds?\b'),
    ('callus_culture', r'callus|cell culture|suspension|protoplast'),
    ('trichome', r'trichome'),
    ('leaf', r'leaf|leaves|leaflet|blade'),
    ('stem_meristem', r'stem|hypocotyl|petiole|pedicel|\bnode|internode|shoot apex|\bsam\b|meristem|axillary|graft'),
    ('seedling_shoot', r'seedlings?|shoot|aerial|whole plant|plantlet|cotyledon'),
]

PROFILES = {
    'rice': {
        'target': 'oryza_sativa',
        'genotype': _rice_genotype,
        'rank': {'japonica_nipponbare': 0, 'japonica_other': 1, 'unknown_sativa': 2, 'indica_or_aus': 3},
        'organs': RICE_ORGANS,
        'deleterious': re.compile(r'male[- ]sterile|\bcms\b|tgms|ptgms|egms|\b\d+s\b|sterile', re.I),
        'extra_mutant': None,
        'extra_other_line': None,
    },
    'tomato': {
        'target': 'solanum_lycopersicum',
        'genotype': _tomato_genotype,
        'rank': {'heinz_1706': 0, 'cultivar_other': 1, 'unknown_lycopersicum': 2, 'cerasiforme': 3,
                 'interspecific_line': 4},
        'organs': TOMATO_ORGANS,
        'deleterious': re.compile(r'male[- ]sterile|\bcms\b|sterile', re.I),
        # tomato ripening/hormone mutants and silencing-pathway genes (sl- prefix); checked on genotype text only
        'extra_mutant': re.compile(r'(?<![a-z])(rin|nor|cnr|flacca|sitiens|notabilis|never ?ripe|green[- ]ripe|'
                                   r'sldcl\d?\w*|slago\d?\w*|slrdr\d?\w*|slrdr|sldrm|slmet1|slcmt)(?![a-z])', re.I),
        'extra_other_line': TOMATO_IL,
    },
}
