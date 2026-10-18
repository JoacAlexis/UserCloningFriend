
import random
import streamlit as st
import html

st.set_page_config(page_title="USER Fusion", layout="wide")

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
            "sequence": "",
            "color": COLORS[0],
            "overlap_preference": "only this seq",
        }
    }


# -------------------------
# Sequence operations
# -------------------------

def add_sequence():
    global color_index
    if not color_index:
        color_index=0
    elif color_index > len(COLOR):
        color_index = 0
    else:
        color_index += 1
    seq_id = st.session_state.next_id
    st.session_state.next_id += 1
    st.session_state.sequence_ids.append(seq_id)

    st.session_state.sequence_data[seq_id] = {
        "gene_name": f"Gene {seq_id}",
        "sequence": "",
        "color": COLORS[color_index],
        "overlap_preference": "only this seq",
    }
    


def remove_sequence(seq_id):
    global color_index
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
    color_index -= 1


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
    primer_length = st.number_input("primer sequence overlap",1,16,8)
    primer_overhang = st.number_input("primer overhang",1,16,8)


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
def MakePrimers(seq1,seq2, overlap):
    pass


for i in range(len(ordered_sequences)-1):
    pass


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
