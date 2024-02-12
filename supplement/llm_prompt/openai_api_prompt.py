"""
This is a sample code exported from OpenAI's API Playground. \
Original results were not obtained through the official API \
but through the Playground (on 2023-04-07). Kept for documentation/reference.
"""
import os
import openai
import json

openai.api_key = ""  # set to own API key

response = openai.Completion.create(
    model="text-davinci-003",
    prompt=(
        # description of the task
        "Given the current state of a list and a prompt, extrapolate as many reactions "
        "of molecules as possible from the prompt and update the list. Every line in "
        "the list contains gene names which act as the reactants and type of reaction "
        "of the left and products on the right, connected with different types of "
        "arrows depending on the reversibility of the reaction. "
        "The available reactions are:\nA dimerizes <--> A2\n"
        "A binds B <--> AB\nAB dissociates to A and B\n"
        "uA is phosphorylated --> pA\npA is dephosphorylated --> uA\n"
        "B phosphorylates uA --> pA\nB dephosphorylates pA --> uA\n"
        "B transcribes a\nB synthesizes A\nA is synthesized\n"
        "B degrades A\nA is degraded\n"
        "Acyt translocates from cytoplasm to nucleus <--> Anuc\n\n\n"
        # one-shot example
        "Examples:\n"
        "current state:\n"
        "[\nA binds B <--> AB\n"
        "AB phosphorylates C --> pC\n"
        "pC is dephosphorylated --> C\n"
        "]\n\n"
        "prompt: The phosphorylated protein C exhibits kinase activity and "
        "phosphorylates its substrate D. Transcription factor D is then translocated "
        "into the nucleaus to transcribe E and F. After E and F are translated into "
        "G and H, respectively, they form a complex (denoted GH) and acquires"
        "phosphatase activity. GH dephosphorylates D, which completes a negative "
        "feedback loop. We also include the dissociation reaction of GH, which yields "
        "free G and H. It is assumed C and D are dephosphorylated at a certain rate.\n"
        "\nnew state:\n"
        "[\nA binds B <--> AB\n"
        "AB phosphorylates C --> pC\n"
        "pC phosphorylates D --> pD\n"
        "pD is translocated to nucleus <--> pDn\n"
        "pDn transcribes E and F\n"
        "E is translated into protein G --> G\n"
        "F is translated into protein H --> H\n"
        "G binds H <--> GH\n"
        "GH dephosphorylates pD --> D\n"
        "GH dissociates to G and H\n"
        "pC is dephosphorylated --> C\n"
        "pD is dephosphorylated --> D\n"
        "]\n\n"
        # input of the task
        "current state:\n"
        "[\nEGF binds EGFR <--> Ra\n"
        "Ra dimerizes <--> R2\n"
        "R2 is phosphorylated --> RP\n"
        "RP is dephosphorylated --> R2\n"
        "]\n\n"
        "prompt: In MedB-1 model IL13 binds to Rec and the resulting species IL13_Rec "
        "triggers the phosphorylation of JAK2. pJAK2 phosphorylates IL13_Rec that "
        "in turn also contributes to the phosphorylation of JAK2. pJAK2 phosphorylates "
        "STAT5 and pSTAT5 activates the transcription of target genes (CD274mRNA, "
        "SOCS3mRNA). The signal strength depends on Rec_i (receptor recycle), on "
        "p_IL13_Rec internalization and degradation and on the negative regulators "
        "SHP1, which de-phosphorylates pJAK2 and pSTAT5, as well as SOCS3, which is "
        "translated from SOCS3mRNA and inhibits the phosphorylation of JAK2. IL13 can "
        "also bind to DecoyR, forming the complex IL13_DecoyR without contributing to "
        "JAK2/STAT5 signaling. In the L1236 model SOCS3mRNA, SOCS3 protein, DecoyR and "
        "the associated reactions are not present.\n\n"
        "new state:\n"
    ),
    temperature=0,
    max_tokens=1024,
    top_p=1,
    frequency_penalty=0,
    presence_penalty=0,
)

if not os.path.exists("out"):
    os.makedirs("out")

# dump response to json
with open("out/response.json", "w") as f:
    json.dump(response, f, indent=4)

# store content
content = response["choices"][0]["message"]["content"]
with open("out/output.txt", "w") as f:
    f.write(content)
