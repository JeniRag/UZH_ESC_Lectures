from Bio import ExPASy, SeqIO

for sid in ['B4F440','O03169','P00850','Q33823','Q95A26',
            'F8RBX8','O21402','P33507','Q35920','Q9T9W0',
            'H6TQD1','P00846','Q2I3G9','Q38PR7','Q9TA24']:
    try:
        handle = ExPASy.get_sprot_raw(sid)
        seq = SeqIO.read(handle,'swiss')
        SeqIO.write(seq, sid+'.genbank','genbank')
        SeqIO.write(seq, sid+'.fasta','fasta')
        print(sid,'sequence length',len(seq))
    except Exception:
        print (sid,'not found')


d={"B4F440": "Neanderthal", "O03169": "Coelacanth", "P00850":"Fruit fly", " Q33823":"Starfish", "Q95A26": " Bornean orangutan",
   "F8RBX8": "Pyhton","O21402": "Ostrich","P33507": "Malaria mosquito", "Q35920": "Atlantic salmon", "Q9T9W0" : "Chimpanze",
   "H6TQD1": "Great fruit-eating bat", "P00846": "Human","Q2I3G9": "Indian elephant", "Q38PR7": "Siberian woolly mammoth", "Q9TA24":"African elephant" }