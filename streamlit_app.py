
import streamlit as st
import html

st.set_page_config(page_title="USER Fusion", layout="wide")

DNA_BASES = "ACGTRYMKBDHVN"
DNA_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVN",
    "TGCAYRKMVHDBN",
)

OVERLAP_OPTIONS = [
    "only this seq",
    "overlapping primer",
    "only next primer",
]

COLORS = [
    "#386CB0",
    "#E07835",
    "#398564",
    "#C94C4C",
    "#8764A8",
    "#287F8E",
    "#B85C83",
]


def normalize_sequence(sequence):
    normalized = "".join(sequence.split()).upper()
    invalid_bases = sorted(set(normalized) - set(DNA_BASES))
    if invalid_bases:
        raise ValueError(
            "Sequence contains unsupported characters: "
            + ", ".join(invalid_bases)
        )
    return normalized


def reverse_complement(sequence):
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def find_overhang(
    seq1,
    seq2,
    preference,
    preferred_length,
    overlap_min=6,
    overlap_max=12,
):
    seq1_window = seq1[-overlap_max * 2:]
    seq2_window = seq2[:overlap_max * 2]

    if preference == "only this seq":
        potential_overhangs = seq1_window
        switch_point = len(seq1_window)
    elif preference == "only next primer":
        potential_overhangs = seq2_window
        switch_point = 0
    else:
        potential_overhangs = seq1_window + seq2_window
        switch_point = len(seq1_window)

    candidates = []
    for start in range(len(potential_overhangs)):
        if potential_overhangs[start] != "A":
            continue
        for length in range(overlap_min, overlap_max + 1):
            end = start + length
            if end > len(potential_overhangs):
                break
            if potential_overhangs[end - 1] != "T":
                continue

            distance_from_junction = min(
                abs(start - switch_point),
                abs(end - switch_point),
            )
            score = (
                -abs(length - preferred_length)
                - 0.5 * distance_from_junction
            )
            candidates.append(
                (score, distance_from_junction, start, potential_overhangs[start:end])
            )

    if not candidates:
        raise ValueError(
            f"No A-to-T overhang between {overlap_min} and "
            f"{overlap_max} bases was found for this junction."
        )

    score, _, _, overhang = max(
        candidates,
        key=lambda candidate: (
            candidate[0],
            -candidate[1],
            -candidate[2],
        ),
    )
    return overhang, score


def MakePrimers(
    seq1,
    seq2,
    overlap_min=6,
    overlap_max=12,
    overhang_preference="overlapping primer",
    preferred_overhang_length=8,
    binding_length=20,
):
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)
    if len(seq1) < binding_length or len(seq2) < binding_length:
        raise ValueError(
            f"Both sequences must be at least {binding_length} bases long "
            "to use the selected primer binding length."
        )

    overhang, score = find_overhang(
        seq1,
        seq2,
        overhang_preference,
        preferred_overhang_length,
        overlap_min,
        overlap_max,
    )
    seq1_binding = seq1[-binding_length:]
    seq2_binding = seq2[:binding_length]

    return {
        "overhang": overhang,
        "score": score,
        "seq1_end": seq1[-(binding_length + overlap_max):],
        "seq2_start": seq2[:binding_length + overlap_max],
        "reverse_primer": reverse_complement(overhang)
        + reverse_complement(seq1_binding),
        "forward_primer": overhang + seq2_binding,
    }


# -------------------------
# Session state
# -------------------------
if "sequence_ids" not in st.session_state:
    st.session_state.sequence_ids = [1]

if "next_id" not in st.session_state:
    st.session_state.next_id = 2

if "sequence_data" not in st.session_state:
    st.session_state.sequence_data = {
        1: {
            "gene_name": "Gene 1",
            "sequence": "ATATATATATATATATAT",
            "color": COLORS[0],
            "overlap_preference": "only this seq",
        }
    }


# -------------------------
# Sequence operations
# -------------------------
if "color_index" not in st.session_state:
    st.session_state.color_index = 1

def add_sequence():
    idx = st.session_state.color_index
    st.session_state.color_index = (idx+1) % len(COLORS) 
    seq_id = st.session_state.next_id
    st.session_state.next_id += 1
    st.session_state.sequence_ids.append(seq_id)

    st.session_state.sequence_data[seq_id] = {
        "gene_name": f"Gene {seq_id}",
        "sequence": "",
        "color": COLORS[idx],
        "overlap_preference": "only this seq",
    }
    
    


def remove_sequence(seq_id):
    idx = st.session_state.color_index
    st.session_state.color_index = (idx+1) % len(COLORS) 
    if len(st.session_state.sequence_ids) <= 1:
        return

    st.session_state.sequence_ids.remove(seq_id)
    st.session_state.sequence_data.pop(seq_id, None)

    # Remove the deleted sequence's widget state
    for key in [
        f"name_{seq_id}",
        f"sequence_{seq_id}",
        f"color_{seq_id}",
        f"overlap_{seq_id}",
    ]:
        st.session_state.pop(key, None)
    


def move_sequence(seq_id, direction):
    ids = st.session_state.sequence_ids
    index = ids.index(seq_id)
    new_index = index + direction

    if 0 <= new_index < len(ids):
        ids[index], ids[new_index] = (
            ids[new_index],
            ids[index],
        )


# -------------------------
# Sidebar styling and inputs
# -------------------------
with st.sidebar:
    st.markdown("## Setup Parameters")
    primer_length = st.number_input(
        "Primer binding length",
        min_value=8,
        max_value=40,
        value=20,
    )
    primer_overhang = st.number_input(
        "Preferred overhang length",
        min_value=6,
        max_value=12,
        value=8,
    )


    st.markdown("## Input sequences")

    st.caption(
        f"{len(st.session_state.sequence_ids)} sequences"
    )

    st.button(
        "+ Add sequence",
        key="add_sequence",
        type="tertiary",
        on_click=add_sequence,
    )

    st.divider()

    for order, seq_id in enumerate(
        st.session_state.sequence_ids, start=1
    ):
        data = st.session_state.sequence_data[seq_id]
        color = data["color"]

        # Color the input fields for this gene
        st.markdown(
            f"""
            <style>
            .st-key-gene_{seq_id} input,
            .st-key-gene_{seq_id} textarea {{
                background-color: {color}18;
            }}

            .st-key-gene_{seq_id} [data-baseweb="select"] > div {{
                background-color: {color}18;
            }}

            .st-key-gene_{seq_id} [data-testid="stExpander"]
            details > summary {{
                background-color: {color}20;
                border-left: 4px solid {color};
                border-radius: 5px;
            }}
            </style>
            """,
            unsafe_allow_html=True,
        )

        with st.container(key=f"gene_{seq_id}"):

            # Fixed-width controls and flexible gene name
            col_up, col_down, col_gene, col_remove = st.columns(
                [0.7, 0.7, 5.6, 0.7],
                gap="small",
            )

            with col_up:
                st.button(
                    "↑",
                    key=f"up_{seq_id}",
                    help="Move gene up",
                    disabled=(order == 1),
                    on_click=move_sequence,
                    args=(seq_id, -1),
                    use_container_width=True,
                )

            with col_down:
                st.button(
                    "↓",
                    key=f"down_{seq_id}",
                    help="Move gene down",
                    disabled=(
                        order == len(st.session_state.sequence_ids)
                    ),
                    on_click=move_sequence,
                    args=(seq_id, 1),
                    use_container_width=True,
                )

            with col_gene:
                # The gene name is now the expander header
                with st.expander(
                    f"{order}. {data['gene_name']}",
                    expanded=False,
                ):
                    data["gene_name"] = st.text_input(
                        "Gene name",
                        value=data["gene_name"],
                        key=f"name_{seq_id}",
                    )

                    data["sequence"] = st.text_area(
                        "DNA sequence",
                        value=data["sequence"],
                        key=f"sequence_{seq_id}",
                        height=120,
                        placeholder="Paste DNA sequence...",
                    )

                    data["color"] = st.color_picker(
                        "Gene color",
                        value=data["color"],
                        key=f"color_{seq_id}",
                    )

                    data["overlap_preference"] = st.selectbox(
                        "Overlap preference",
                        OVERLAP_OPTIONS,
                        index=OVERLAP_OPTIONS.index(
                            data["overlap_preference"]
                        ),
                        key=f"overlap_{seq_id}",
                    )

            with col_remove:
                st.button(
                    "×",
                    key=f"remove_{seq_id}",
                    help="Remove gene",
                    type="tertiary",
                    disabled=(
                        len(st.session_state.sequence_ids) <= 1
                    ),
                    on_click=remove_sequence,
                    args=(seq_id,),
                    use_container_width=True,
                )

    # Rebuild data in the current assembly order
    ordered_sequences = [
        {
            "id": seq_id,
            "order": order,
            **st.session_state.sequence_data[seq_id],
        }
        for order, seq_id in enumerate(
            st.session_state.sequence_ids, start=1
        )
    ]

   


# -------------------------
# Make primers
# -------------------------
st.subheader("Junction primers")

if len(ordered_sequences) < 2:
    st.info("Add at least two sequences to design junction primers.")

for index in range(len(ordered_sequences) - 1):
    left = ordered_sequences[index]
    right = ordered_sequences[index + 1]
    junction_name = f"{left['gene_name']} → {right['gene_name']}"

    with st.expander(junction_name, expanded=True):
        left_sequence = left["sequence"]
        right_sequence = right["sequence"]
        if not left_sequence.strip() or not right_sequence.strip():
            st.info("Enter both sequences to design primers for this junction.")
            continue

        try:
            primers = MakePrimers(
                left_sequence,
                right_sequence,
                overhang_preference=left["overlap_preference"],
                preferred_overhang_length=primer_overhang,
                binding_length=primer_length,
            )
        except ValueError as error:
            st.warning(str(error))
            continue

        st.caption(
            f"Selected A-to-T overhang ({len(primers['overhang'])} bp): "
            f"{primers['overhang']} · score {primers['score']:.1f}"
        )
        left_column, right_column = st.columns(2)
        with left_column:
            st.markdown(f"**{left['gene_name']} end**")
            st.code(primers["seq1_end"], language=None)
            st.markdown("**Reverse primer (5′→3′)**")
            st.code(primers["reverse_primer"], language=None)
        with right_column:
            st.markdown(f"**{right['gene_name']} start**")
            st.code(primers["seq2_start"], language=None)
            st.markdown("**Forward primer (5′→3′)**")
            st.code(primers["forward_primer"], language=None)


# -------------------------
# Main page
# -------------------------
st.title("USER Fusion Designer")

LINE_WIDTH = 80

st.subheader("Assemble Order")
assembly_html = []
line_position = 0

for seq_id in st.session_state.sequence_ids:
    data = st.session_state.sequence_data[seq_id]

    gene_name = data["gene_name"]
    #sequence = "".join(data["sequence"].split())
    color = data["color"]

    #if not sequence:
    #    continue
    # Split each gene as needed to keep the entire assembly
    # wrapped at exactly LINE_WIDTH characters.
    #position = 0

    assembly_html.append(
        f'<span title="{html.escape(gene_name)}" '
        f'style="background-color:{html.escape(color)}; '
        f'color:#111; padding:1px 0;">'
        f'{html.escape(gene_name)}</span>'
        )

st.markdown(
f"""
<div style="
     font-family: monospace;
     font-size: 15px;
     line-height: 1.8;
     overflow-x: auto;
     overflow-wrap: normal;
     white-space: normal;
     max-width: 80ch;
 ">
{''.join(assembly_html)}
</div>
""",
unsafe_allow_html=True,
)



st.subheader("Full assembly")

assembly_html = []
line_position = 0

for seq_id in st.session_state.sequence_ids:
    data = st.session_state.sequence_data[seq_id]

    gene_name = data["gene_name"]
    sequence = "".join(data["sequence"].split())
    color = data["color"]

    if not sequence:
        continue
    # Split each gene as needed to keep the entire assembly
    # wrapped at exactly LINE_WIDTH characters.
    position = 0

    while position < len(sequence):
        available = LINE_WIDTH - line_position
        chunk = sequence[position:position + available]

        assembly_html.append(
            f'<span title="{html.escape(gene_name)}" '
            f'style="background-color:{html.escape(color)}; '
            f'color:#111; padding:1px 0;">'
            f'{html.escape(chunk)}</span>'
        )

        position += len(chunk)
        line_position += len(chunk)

        if line_position == LINE_WIDTH:
            assembly_html.append("<br>")
            line_position = 0
st.markdown(
f"""
<div style="
     font-family: monospace;
     font-size: 15px;
     line-height: 1.8;
     overflow-x: auto;
     overflow-wrap: normal;
     white-space: normal;
     max-width: 80ch;
 ">
{''.join(assembly_html)}
</div>
""",
unsafe_allow_html=True,
)



st.write("Current assembly order")

for item in ordered_sequences:
    st.markdown(
        f"<span style='display:inline-block;"
        f"background:{item['color']};color:white;"
        f"padding:5px 10px;border-radius:5px;"
        f"margin:3px'>"
        f"{item['gene_name']}"
        f"</span>",
        unsafe_allow_html=True,
    )

with st.expander("Inspect sequence data"):
    st.json(ordered_sequences)
