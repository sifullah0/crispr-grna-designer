import re
import RNA
import heapq
import azimuth.model_comparison
import pandas as pd
import numpy as np
import warnings
import sys
import os
from contextlib import contextmanager
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import mysql.connector
from mysql.connector import Error
import datetime

@contextmanager
def suppress_stdout():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

warnings.filterwarnings("ignore")

def get_db_connection():
    try:
        connection = mysql.connector.connect(
            host='localhost',  # Adjust as needed
            database='grna_db',  # Create this database manually
            user='root',  # Adjust as needed
            password='root'  # Adjust as needed
        )
        return connection
    except Error as e:
        messagebox.showerror("Database Error", f"Failed to connect to database: {str(e)}")
        return None

def initialize_db():
    connection = get_db_connection()
    if connection:
        cursor = connection.cursor()
        try:
            # Create sequences table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS sequences (
                    id VARCHAR(50) PRIMARY KEY,
                    sequence TEXT NOT NULL,
                    timestamp DATETIME NOT NULL
                )
            """)
            # Create results table with ON DELETE CASCADE
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS results (
                    id INT AUTO_INCREMENT PRIMARY KEY,
                    sequence_id VARCHAR(50) NOT NULL,
                    grna_sequence VARCHAR(20) NOT NULL,
                    start INT NOT NULL,
                    end INT NOT NULL,
                    azimuth_score FLOAT NOT NULL,
                    accessibility_score FLOAT NOT NULL,
                    total_score FLOAT NOT NULL,
                    mfe FLOAT NOT NULL,
                    structure VARCHAR(255) NOT NULL,
                    pam_sequence VARCHAR(3) NOT NULL,
                    gc_content FLOAT NOT NULL,
                    FOREIGN KEY (sequence_id) REFERENCES sequences(id) ON DELETE CASCADE
                )
            """)
            connection.commit()
        except Error as e:
            messagebox.showerror("Database Error", f"Failed to initialize database: {str(e)}")
        finally:
            cursor.close()
            connection.close()

# Call initialize_db at startup
initialize_db()

def find_forward_grna_candidates(target_sequence, pam_pattern=r'[ACGT]GG'):
    seq = target_sequence.upper()
    candidates = []
    for i in range(len(seq) - 22):
        window = seq[i : i + 23]
        spacer, pam = window[:20], window[20:]
        if re.fullmatch(pam_pattern, pam):
            candidates.append({
                "sequence": spacer,
                "pam": pam,
                "start": i + 1,
                "end": i + 23
            })
    return candidates

def calculate_gc_content(gRNA_seq):
    gc_count = gRNA_seq.count('G') + gRNA_seq.count('C')
    return gc_count / len(gRNA_seq)

def filter_grna(gRNA_seq, gc_min=0.4, gc_max=0.6):
    gc_content = calculate_gc_content(gRNA_seq)
    if not (gc_min <= gc_content <= gc_max):
        return False
    if 'TTTT' in gRNA_seq:
        return False
    return True

def gRNA_score(gRNA_seq):
    fc = RNA.fold_compound(gRNA_seq)
    structure, mfe = fc.mfe()
    if mfe < -15:
        return ('mfe', mfe, structure)
    else:
        fc.exp_params_rescale(mfe)
        fc.pf()
        plist = fc.plist_from_probs(1e-6)
        bpp = {(ep.i, ep.j): ep.p for ep in plist}
        unpaired_probs = []
        for i in range(1, min(13, len(gRNA_seq) + 1)):
            paired = sum(p for (a, b), p in bpp.items() if a == i or b == i)
            unpaired_probs.append(1.0 - paired)
        avg_unpaired = sum(unpaired_probs) / len(unpaired_probs)
        return (avg_unpaired, structure, mfe)

def get_azimuth_score(gRNA_seq, target_sequence, start):
    seq = target_sequence.upper()
    start_idx = max(0, start - 5)
    end_idx = min(len(seq), start + 25)
    sequence_30mer = seq[start_idx:end_idx].ljust(30, 'N')
    data = pd.DataFrame({'sequence': [sequence_30mer]})
    try:
        with suppress_stdout():
            scores = azimuth.model_comparison.predict(data['sequence'].to_numpy())
        return scores[0]
    except:
        return 0.0

class gRNAApp:
    def __init__(self, root):
        self.root = root
        self.root.title("gRNA Candidate Finder")
        self.root.geometry("1200x700")
        
        # Configure styles first
        style = ttk.Style()
        style.configure('Bold.TLabelframe.Label', font=('Arial', 12, 'bold'))
        style.configure('Bold.TButton', font=('Arial', 10, 'bold'))
        style.configure("Bold.Treeview.Heading", font=("Arial", 10, "bold"))

        # Main container with padding
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Input frame with better organization
        input_frame = ttk.LabelFrame(main_frame, text="Input Parameters", padding="10", style='Bold.TLabelframe')
        input_frame.pack(fill=tk.X, pady=(0, 10))
        
        # DNA Sequence input - wider and with scrollbar
        ttk.Label(input_frame, text="DNA Sequence:", font=('Arial', 10, 'bold')).grid(row=0, column=0, sticky=tk.W, pady=(0, 5))
        
        seq_frame = ttk.Frame(input_frame)
        seq_frame.grid(row=1, column=0, columnspan=4, sticky=(tk.W, tk.E), pady=(0, 10))
        
        seq_scrollbar_y = ttk.Scrollbar(seq_frame, orient=tk.VERTICAL)
        
        self.seq_entry = tk.Text(seq_frame, height=5, wrap=tk.CHAR, 
                                  yscrollcommand=seq_scrollbar_y.set,
                                  font=('Courier', 10))
        
        seq_scrollbar_y.config(command=self.seq_entry.yview)
        
        self.seq_entry.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        seq_scrollbar_y.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        seq_frame.columnconfigure(0, weight=1)
        seq_frame.rowconfigure(0, weight=1)
        
        # Parameters in a grid layout
        params_frame = ttk.Frame(input_frame)
        params_frame.grid(row=2, column=0, columnspan=4, sticky=(tk.W, tk.E))
        
        ttk.Label(params_frame, text="PAM Pattern:", font=('Arial', 9, 'bold')).grid(row=0, column=0, sticky=tk.W, padx=(0, 10))
        self.pam_entry = ttk.Entry(params_frame, width=20)
        self.pam_entry.insert(0, "[ACGT]GG")
        self.pam_entry.grid(row=0, column=1, sticky=tk.W, padx=(0, 30))
        
        ttk.Label(params_frame, text="GC Content Min:", font=('Arial', 9, 'bold')).grid(row=0, column=2, sticky=tk.W, padx=(0, 10))
        self.gc_min = ttk.Entry(params_frame, width=10)
        self.gc_min.insert(0, "0.4")
        self.gc_min.grid(row=0, column=3, sticky=tk.W, padx=(0, 30))
        
        ttk.Label(params_frame, text="GC Content Max:", font=('Arial', 9, 'bold')).grid(row=0, column=4, sticky=tk.W, padx=(0, 10))
        self.gc_max = ttk.Entry(params_frame, width=10)
        self.gc_max.insert(0, "0.6")
        self.gc_max.grid(row=0, column=5, sticky=tk.W)
        
        # Run button - larger and more prominent
        ttk.Button(input_frame, text="Find gRNA Candidates", padding="4",
                   command=self.run_analysis, 
                   style='Bold.TButton').grid(row=3, column=0, columnspan=4, pady=(20, 10))
        
        input_frame.columnconfigure(0, weight=1)
        
        # Results frame with scrollbars
        results_frame = ttk.LabelFrame(main_frame, text="Results", padding="10", style='Bold.TLabelframe')
        results_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create frame for treeview and scrollbars
        tree_frame = ttk.Frame(results_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))
        
        # Scrollbars
        vsb = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL)
        
        # Treeview with both scrollbars
        self.tree = ttk.Treeview(tree_frame, 
                                 columns=("Sequence", "Start", "End", "Azimuth", "Accessibility", "Total", "MFE", "Structure"), 
                                 show="headings",
                                 yscrollcommand=vsb.set,
                                 height=12,
                                 style="Bold.Treeview")
        
        vsb.config(command=self.tree.yview)
        
        # Configure columns with appropriate widths
        self.tree.heading("Sequence", text="Sequence")
        self.tree.heading("Start", text="Start")
        self.tree.heading("End", text="End")
        self.tree.heading("Azimuth", text="Azimuth Score")
        self.tree.heading("Accessibility", text="Accessibility")
        self.tree.heading("Total", text="Total Score")
        self.tree.heading("MFE", text="MFE")
        self.tree.heading("Structure", text="Structure")
        
        self.tree.column("Sequence", width=180, anchor=tk.CENTER)
        self.tree.column("Start", width=60, anchor=tk.CENTER)
        self.tree.column("End", width=60, anchor=tk.CENTER)
        self.tree.column("Azimuth", width=100, anchor=tk.CENTER)
        self.tree.column("Accessibility", width=100, anchor=tk.CENTER)
        self.tree.column("Total", width=100, anchor=tk.CENTER)
        self.tree.column("MFE", width=80, anchor=tk.CENTER)
        self.tree.column("Structure", width=200, anchor=tk.CENTER)
        
        # Grid layout for tree and scrollbars
        self.tree.grid(row=0, column=0, sticky=(tk.N, tk.S, tk.E, tk.W))
        vsb.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.rowconfigure(0, weight=1)
        
        # Button frame for Save CSV, Save DB, Visualize, and Manage DB buttons
        button_frame = ttk.Frame(results_frame)
        button_frame.pack(side=tk.RIGHT, pady=(8, 0))
        
        ttk.Button(button_frame, text="Visualize Results", padding="4", 
                   command=self.visualize_results, style='Bold.TButton').pack(side=tk.LEFT, padx=(0, 5))
        
        ttk.Button(button_frame, text="Save Results to CSV", padding="4", 
                   command=self.save_results_csv, style='Bold.TButton').pack(side=tk.LEFT, padx=(0, 5))
        
        ttk.Button(button_frame, text="Save to Database", padding="4", 
                   command=self.save_to_database, style='Bold.TButton').pack(side=tk.LEFT, padx=(0, 5))
        
        ttk.Button(button_frame, text="Manage Database", padding="4", 
                   command=self.manage_database, style='Bold.TButton').pack(side=tk.LEFT)

        # Store scored_list for later use
        self.scored_list = []

    def run_analysis(self):
        self.tree.delete(*self.tree.get_children())
        self.root.update()
        
        sequence = self.seq_entry.get("1.0", tk.END).strip()
        
        if not sequence:
            messagebox.showerror("Error", "Please enter a DNA sequence")
            return
        
        pam_pattern = self.pam_entry.get()
        try:
            gc_min = float(self.gc_min.get())
            gc_max = float(self.gc_max.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid GC content values")
            return
        
        try:
            candidates = find_forward_grna_candidates(sequence, pam_pattern)
            filtered_candidates = [cand for cand in candidates if filter_grna(cand["sequence"], gc_min, gc_max)]

            self.root.update()
            
            scored_list = []
            heap = []
            for cand in filtered_candidates:
                seq = cand["sequence"]
                azimuth_score = get_azimuth_score(seq, sequence, cand["start"])
                result = gRNA_score(seq)
                if result[0] == 'mfe':
                    _, mfe_value, struct = result
                    accessibility_score = mfe_value
                else:
                    accessibility_score, struct, mfe = result
                total_score = 0.6 * azimuth_score + 0.4 * accessibility_score
                gc_content = calculate_gc_content(seq)
                heapq.heappush(heap, (-total_score, seq, cand["start"], cand["end"], struct, accessibility_score, azimuth_score, mfe, cand["pam"], gc_content))
            
            while heap:
                neg_score, seq, start, end, struct, accessibility_score, azimuth_score, mfe, pam, gc_content = heapq.heappop(heap)
                scored_list.append({
                    "sequence": seq,
                    "start": start,
                    "end": end,
                    "accessibility_score": accessibility_score,
                    "azimuth_score": azimuth_score,
                    "total_score": -neg_score,
                    "mfe": mfe,
                    "structure": struct,
                    "pam": pam,
                    "gc_content": gc_content
                })
            
            self.scored_list = scored_list
            
            for entry in scored_list:
                self.tree.insert("", tk.END, values=(
                    entry["sequence"],
                    entry["start"],
                    entry["end"],
                    f"{entry['azimuth_score']:.3f}",
                    f"{entry['accessibility_score']:.3f}",
                    f"{entry['total_score']:.3f}",
                    f"{entry['mfe']:.2f}",
                    entry["structure"]
                ))
            
            messagebox.showinfo("Success", f"Found {len(scored_list)} gRNA candidates")
        
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")
    
    def save_to_database(self):
        if not self.scored_list:
            messagebox.showwarning("Warning", "No results to save")
            return
        
        sequence = self.seq_entry.get("1.0", tk.END).strip()
        
        # Prompt for sequence ID
        seq_id = tk.simpledialog.askstring("Sequence ID", "Enter Sequence ID (alphanumeric):")
        if not seq_id:
            return
        
        connection = get_db_connection()
        if connection:
            cursor = connection.cursor()
            try:
                # Insert sequence with timestamp
                timestamp = datetime.datetime.now()
                cursor.execute("""
                    INSERT INTO sequences (id, sequence, timestamp) 
                    VALUES (%s, %s, %s) 
                    ON DUPLICATE KEY UPDATE sequence = %s, timestamp = %s
                """, (seq_id, sequence, timestamp, sequence, timestamp))
                
                # Insert results
                for entry in self.scored_list:
                    cursor.execute("""
                        INSERT INTO results 
                        (sequence_id, grna_sequence, start, end, azimuth_score, accessibility_score, total_score, mfe, structure, pam_sequence, gc_content)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    """, (seq_id, entry["sequence"], entry["start"], entry["end"], entry["azimuth_score"],
                          entry["accessibility_score"], entry["total_score"], entry["mfe"], entry["structure"],
                          entry["pam"], entry["gc_content"]))
                
                connection.commit()
                messagebox.showinfo("Success", "Results saved to database")
            except Error as e:
                messagebox.showerror("Database Error", f"Failed to save to database: {str(e)}")
            finally:
                cursor.close()
                connection.close()

    def manage_database(self):
        manage_window = tk.Toplevel(self.root)
        manage_window.title("Manage Database")
        manage_window.geometry("900x650")
        manage_window.transient(self.root)
        manage_window.grab_set()

        # ---- Styles ----------------------------------------------------
        style = ttk.Style()
        style.configure('Bold.TNotebook.Tab', font=('Arial', 12, 'bold'))
        style.configure('Query.TButton', font=('Arial', 11, 'bold'), padding=6)

        # ---- Notebook --------------------------------------------------
        notebook = ttk.Notebook(manage_window, style='Bold.TNotebook')
        notebook.pack(fill=tk.BOTH, expand=True, padx=15, pady=15)

        # ------------------- SEARCH QUERIES TAB -------------------------
        search_frame = ttk.Frame(notebook, padding=15)
        notebook.add(search_frame, text="Search Queries")

        search_frame.pack_propagate(False)
        search_frame.grid_rowconfigure(0, weight=1)
        search_frame.grid_columnconfigure(0, weight=1)

        inner = ttk.Frame(search_frame)
        inner.grid(row=0, column=0, sticky="nsew")

        queries = [
            ("All gRNAs for a given sequence",                self.query_all_grnas_for_sequence),
            ("Top-N best gRNAs (by total score)",              self.query_top_n_grnas),
            ("Top-N best gRNAs by score threshold",            self.query_top_n_by_threshold),   # NEW
            ("gRNAs with a specific PAM",                      self.query_grnas_by_pam),
            ("gRNAs in a GC-range",                            self.query_grnas_by_gc_range),
            ("Search by partial DNA",                          self.query_by_partial_dna),
            ("Count gRNAs per sequence",                       self.query_count_grnas_per_sequence),
            ("Export a full analysis for a gene",              self.query_full_analysis_for_gene),
        ]

        for i, (label, func) in enumerate(queries):
            btn = ttk.Button(inner, text=label, command=func,
                             style='Query.TButton', width=50)
            btn.grid(row=i, column=0, pady=7, padx=10, sticky="w")

        # Add vertical scrollbar if needed
        canvas = tk.Canvas(search_frame)
        vbar = ttk.Scrollbar(search_frame, orient="vertical", command=canvas.yview)
        scroll_inner = ttk.Frame(canvas)
        scroll_inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.configure(yscrollcommand=vbar.set)
        canvas.create_window((0,0), window=scroll_inner, anchor="nw")
        canvas.grid(row=0, column=0, sticky="nsew")
        vbar.grid(row=0, column=1, sticky="ns")

        # Move buttons to scroll_inner
        for widget in inner.winfo_children():
            widget.destroy()
        for i, (label, func) in enumerate(queries):
            btn = ttk.Button(scroll_inner, text=label, command=func,
                             style='Query.TButton', width=50)
            btn.grid(row=i, column=0, pady=7, padx=10, sticky="w")

        # ------------------- DELETE OPERATIONS TAB ----------------------
        delete_frame = ttk.Frame(notebook, padding=15)
        notebook.add(delete_frame, text="Delete Operations")

        delete_frame.pack_propagate(False)
        delete_frame.grid_rowconfigure(0, weight=1)
        delete_frame.grid_columnconfigure(0, weight=1)

        del_inner = ttk.Frame(delete_frame)
        del_inner.grid(row=0, column=0, sticky="nsew")

        delete_ops = [
            ("Delete individual gRNA", self.delete_individual_grna),
            ("Delete all for sequence", self.delete_all_for_sequence),
        ]

        for i, (label, func) in enumerate(delete_ops):
            btn = ttk.Button(del_inner, text=label, command=func,
                             style='Query.TButton', width=50)
            btn.grid(row=i, column=0, pady=7, padx=10, sticky="w")
    
    def execute_query(self, sql, params=None, fetch=True):
        connection = get_db_connection()
        if connection:
            cursor = connection.cursor(dictionary=True)
            try:
                cursor.execute(sql, params or ())
                if fetch:
                    results = cursor.fetchall()
                    return results
                else:
                    connection.commit()
            except Error as e:
                messagebox.showerror("Database Error", f"Query failed: {str(e)}")
            finally:
                cursor.close()
                connection.close()
        return []

    def show_query_results(self, results, title):
        if not results:
            messagebox.showinfo("No Results", "No results found")
            return
        
        result_window = tk.Toplevel(self.root)
        result_window.title(title)
        result_window.geometry("1000x600")
        
        tree_frame = ttk.Frame(result_window)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        vsb = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL)
        tree = ttk.Treeview(tree_frame, columns=list(results[0].keys()), show="headings", yscrollcommand=vsb.set)
        vsb.config(command=tree.yview)
        
        for col in tree["columns"]:
            tree.heading(col, text=col.upper())
            tree.column(col, width=150, anchor=tk.CENTER)
        
        for row in results:
            tree.insert("", tk.END, values=tuple(row.values()))
        
        tree.grid(row=0, column=0, sticky=(tk.N, tk.S, tk.E, tk.W))
        vsb.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.rowconfigure(0, weight=1)

    def query_all_grnas_for_sequence(self):
        seq_id = tk.simpledialog.askstring("Query", "Enter sequence ID:")
        if seq_id:
            sql = "SELECT * FROM results WHERE sequence_id = %s"
            results = self.execute_query(sql, (seq_id,))
            self.show_query_results(results, "All gRNAs for Sequence")

    def query_top_n_grnas(self):
        seq_id = tk.simpledialog.askstring("Query", "Enter sequence ID:")
        n = tk.simpledialog.askinteger("Query", "Enter N (top N):", minvalue=1)
        if seq_id and n:
            sql = """
                SELECT id AS grna_id, grna_sequence, total_score 
                FROM results 
                WHERE sequence_id = %s 
                ORDER BY total_score DESC 
                LIMIT %s
            """
            results = self.execute_query(sql, (seq_id, n))
            self.show_query_results(results, f"Top {n} gRNAs for Sequence")

    def query_grnas_by_pam(self):
        pam = tk.simpledialog.askstring("Query", "Enter PAM sequence (e.g., NAG):")
        if pam:
            sql = "SELECT * FROM results WHERE pam_sequence = %s"
            results = self.execute_query(sql, (pam,))
            self.show_query_results(results, f"gRNAs with PAM {pam}")

    def query_grnas_by_gc_range(self):
        min_gc = tk.simpledialog.askfloat("Query", "Enter min GC content:", minvalue=0.0, maxvalue=1.0)
        max_gc = tk.simpledialog.askfloat("Query", "Enter max GC content:", minvalue=0.0, maxvalue=1.0)
        if min_gc is not None and max_gc is not None:
            sql = "SELECT * FROM results WHERE gc_content BETWEEN %s AND %s"
            results = self.execute_query(sql, (min_gc, max_gc))
            self.show_query_results(results, f"gRNAs in GC range {min_gc}-{max_gc}")

    def query_by_partial_dna(self):
        motif = tk.simpledialog.askstring("Query", "Enter partial DNA motif:")
        if motif:
            sql = """
                SELECT s.id, s.sequence, r.grna_sequence 
                FROM sequences s 
                LEFT JOIN results r ON s.id = r.sequence_id 
                WHERE s.sequence LIKE %s
            """
            results = self.execute_query(sql, (f"%{motif}%",))
            self.show_query_results(results, "Search by Partial DNA")

    def query_count_grnas_per_sequence(self):
        sql = """
            SELECT s.id, COUNT(r.id) AS num_grnas 
            FROM sequences s 
            LEFT JOIN results r ON s.id = r.sequence_id 
            GROUP BY s.id
        """
        results = self.execute_query(sql)
        self.show_query_results(results, "Count gRNAs per Sequence")

    def query_top_n_by_threshold(self):
        seq_id = tk.simpledialog.askstring(
            "Sequence ID", "Enter sequence ID:")
        if not seq_id:
            return

        n = tk.simpledialog.askinteger(
            "Top-N", "How many top gRNAs do you want?", minvalue=1, maxvalue=1000)
        if n is None:
            return

        threshold = tk.simpledialog.askfloat(
            "Score Threshold",
            "Minimum total_score (0–1, e.g. 0.7):",
            minvalue=0.0, maxvalue=1.0)
        if threshold is None:
            return

        sql = """
            SELECT id AS grna_id, grna_sequence, total_score, azimuth_score,
                   accessibility_score, mfe, pam_sequence, gc_content
            FROM results
            WHERE sequence_id = %s AND total_score >= %s
            ORDER BY total_score DESC
            LIMIT %s
        """
        results = self.execute_query(sql, (seq_id, threshold, n))
        title = f"Top {n} gRNAs for '{seq_id}' (score ≥ {threshold:.3f})"
        self.show_query_results(results, title)

    def query_full_analysis_for_gene(self):
        gene_id = tk.simpledialog.askstring("Query", "Enter gene/sequence ID (e.g., BRCA1-exon2):")
        if gene_id:
            sql = """
                SELECT s.id AS name, r.* 
                FROM sequences s 
                JOIN results r ON s.id = r.sequence_id 
                WHERE s.id = %s
            """
            results = self.execute_query(sql, (gene_id,))
            self.show_query_results(results, f"Full Analysis for {gene_id}")

    def delete_individual_grna(self):
        grna_id = tk.simpledialog.askinteger("Delete", "Enter gRNA ID to delete:")
        if grna_id:
            sql = "DELETE FROM results WHERE id = %s"
            self.execute_query(sql, (grna_id,), fetch=False)
            messagebox.showinfo("Success", f"gRNA {grna_id} deleted")

    def delete_all_for_sequence(self):
        seq_id = tk.simpledialog.askstring("Delete", "Enter sequence ID to delete all gRNAs:")
        if seq_id:
            sql = "DELETE FROM sequences WHERE id = %s"  # Cascade will delete results
            self.execute_query(sql, (seq_id,), fetch=False)
            messagebox.showinfo("Success", f"All data for sequence {seq_id} deleted")

    def visualize_results(self):
        if not self.tree.get_children():
            messagebox.showwarning("Warning", "No results to visualize")
            return
        
        # Extract data from tree
        sequences = []
        azimuth_scores = []
        accessibility_scores = []
        mfe_scores = []
        total_scores = []
        
        for item in self.tree.get_children():
            values = self.tree.item(item)["values"]
            sequences.append(values[0][:10] + "...")  # Truncate sequence for x-axis
            azimuth_scores.append(float(values[3]))
            accessibility_scores.append(float(values[4]))
            total_scores.append(float(values[5]))
            mfe_scores.append(float(values[6]))
        
        total_candidates = len(sequences)
        
        # If more than 15 candidates, ask user how many to display
        if total_candidates > 15:
            choice_window = tk.Toplevel(self.root)
            choice_window.title("Choose Visualization Range")
            choice_window.geometry("400x200")
            choice_window.transient(self.root)
            choice_window.grab_set()
            
            ttk.Label(choice_window, text=f"Found {total_candidates} candidates.", 
                     font=('Arial', 11, 'bold')).pack(pady=(20, 10))
            ttk.Label(choice_window, text="How many would you like to visualize?", 
                     font=('Arial', 10)).pack(pady=(0, 20))
            
            num_to_show = tk.IntVar(value=15)
            
            button_frame = ttk.Frame(choice_window)
            button_frame.pack(pady=10)
            
            def show_selection():
                choice_window.destroy()
                self._create_visualization(sequences[:num_to_show.get()], 
                                          azimuth_scores[:num_to_show.get()],
                                          accessibility_scores[:num_to_show.get()],
                                          mfe_scores[:num_to_show.get()],
                                          total_scores[:num_to_show.get()],
                                          total_candidates)
            
            ttk.Button(button_frame, text="Top 15", padding="5",
                      command=lambda: [num_to_show.set(15), show_selection()],
                      style='Bold.TButton').pack(side=tk.LEFT, padx=5)
            
            ttk.Button(button_frame, text="Top 30", padding="5",
                      command=lambda: [num_to_show.set(min(30, total_candidates)), show_selection()],
                      style='Bold.TButton').pack(side=tk.LEFT, padx=5)
            
            ttk.Button(button_frame, text="All", padding="5",
                      command=lambda: [num_to_show.set(total_candidates), show_selection()],
                      style='Bold.TButton').pack(side=tk.LEFT, padx=5)
            
            # Custom number entry
            custom_frame = ttk.Frame(choice_window)
            custom_frame.pack(pady=10)
            ttk.Label(custom_frame, text="Custom number:").pack(side=tk.LEFT, padx=5)
            custom_entry = ttk.Entry(custom_frame, width=10)
            custom_entry.pack(side=tk.LEFT, padx=5)
            
            def show_custom():
                try:
                    num = int(custom_entry.get())
                    if 1 <= num <= total_candidates:
                        num_to_show.set(num)
                        show_selection()
                    else:
                        messagebox.showerror("Error", f"Please enter a number between 1 and {total_candidates}")
                except ValueError:
                    messagebox.showerror("Error", "Please enter a valid number")
            
            ttk.Button(custom_frame, text="Show", padding="5",
                      command=show_custom).pack(side=tk.LEFT, padx=5)
            
            choice_window.wait_window()
        else:
            self._create_visualization(sequences, azimuth_scores, accessibility_scores, 
                                      mfe_scores, total_scores, total_candidates)
    
    def _create_visualization(self, sequences, azimuth_scores, accessibility_scores, 
                             mfe_scores, total_scores, total_candidates):
        # Create new window for visualization
        viz_window = tk.Toplevel(self.root)
        viz_window.title(f"gRNA Candidates Visualization (Showing {len(sequences)} of {total_candidates})")
        viz_window.geometry("1400x800")
        
        
        
        # Create main container with both scrollbars
        main_container = ttk.Frame(viz_window)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=(5, 10))
        
        # Create canvas with both scrollbars
        canvas_widget = tk.Canvas(main_container, highlightthickness=0)
        v_scrollbar = ttk.Scrollbar(main_container, orient=tk.VERTICAL, command=canvas_widget.yview)
        h_scrollbar = ttk.Scrollbar(main_container, orient=tk.HORIZONTAL, command=canvas_widget.xview)
        scrollable_frame = ttk.Frame(canvas_widget)
        
        canvas_widget.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Pack scrollbars and canvas
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas_widget.pack( fill=tk.BOTH, expand=True, pady=(10, 0))
        
        # Add save plot button at the top - outside scrollable area
        button_frame = ttk.Frame(viz_window)
        button_frame.pack(side=tk.BOTTOM, pady=(5, 20))
        ttk.Button(button_frame, text="Save Plot as PNG", padding="4",
                   command=lambda: self.save_plot(fig), 
                   style='Bold.TButton').pack()
        
        # Adjust figure size based on number of candidates
        fig_width = max(12, len(sequences) * 0.6)
        fig_height = 8  # Increased height to accommodate legend below
        fig = Figure(figsize=(fig_width, fig_height), dpi=100)
        ax = fig.add_subplot(111)
        
        # Set up bar positions
        x = np.arange(len(sequences))
        width = 0.2
        
        # Create grouped bars
        bars1 = ax.bar(x - 1.5*width, azimuth_scores, width, label='Azimuth Score', color='#3498db', alpha=0.8)
        bars2 = ax.bar(x - 0.5*width, accessibility_scores, width, label='Accessibility Score', color='#2ecc71', alpha=0.8)
        bars3 = ax.bar(x + 0.5*width, mfe_scores, width, label='MFE', color='#e74c3c', alpha=0.8)
        bars4 = ax.bar(x + 1.5*width, total_scores, width, label='Total Score', color='#f39c12', alpha=0.8)
        
        # Customize plot
        ax.set_xlabel('gRNA Candidates', fontsize=12, fontweight='bold', labelpad=10)
        ax.set_ylabel('Score Values', fontsize=12, fontweight='bold')
        title = f'gRNA Candidate Scores Comparison (Top {len(sequences)} of {total_candidates})' if len(sequences) < total_candidates else f'gRNA Candidate Scores Comparison (All {total_candidates} Candidates)'
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        ax.set_xticks(x)
        ax.set_xticklabels(sequences, rotation=45, ha='right', fontsize=max(6, 10 - len(sequences)//10))
        
        # Set y-axis ticks at intervals of 0.2
        all_scores = np.concatenate([azimuth_scores, accessibility_scores, mfe_scores, total_scores])
        min_score = np.floor(min(all_scores) / 0.4) * 0.4  # Round down to nearest 0.4
        max_score = np.ceil(max(all_scores) / 0.4) * 0.4   # Round up to nearest 0.4
        y_ticks = np.arange(min_score, max_score + 0.4, 0.4)  # Ticks at 0.4 intervals
        ax.set_yticks(y_ticks)
        ax.set_ylim(min_score - 0.1, max_score + 0.1)  # Add small margin
        
        # Place legend below the plot with more space from x-axis label
        ax.legend(loc='lower left', bbox_to_anchor=(0.86, 1.01),
                 ncol=1, fontsize=8, frameon=True, shadow=True)
        
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.axhline(y=0, color='black', linewidth=0.8)
        
        # Adjust layout to prevent label cutoff and accommodate legend with more bottom space
        fig.tight_layout(rect=[0, 0, 1, 1])
        
        # Embed plot in scrollable frame
        canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()
        
        # Configure scrollable region
        canvas_widget.create_window((0, 0), window=scrollable_frame, anchor=tk.NW)
        scrollable_frame.update_idletasks()
        canvas_widget.config(scrollregion=canvas_widget.bbox("all"))
    
    def save_results_csv(self):
        if not self.tree.get_children():
            messagebox.showwarning("Warning", "No results to save")
            return
        
        import csv
        try:
            with open("grna_candidates.csv", "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Sequence", "Start", "End", "Azimuth Score", "Accessibility Score", "Total Score", "MFE", "Structure"])
                for item in self.tree.get_children():
                    writer.writerow(self.tree.item(item)["values"])
            messagebox.showinfo("Success", "Results saved to grna_candidates.csv")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save results: {str(e)}")
    
    def save_plot(self, fig):
        try:
            fig.savefig("grna_visualization.png")
            messagebox.showinfo("Success", "Plot saved as grna_visualization.png")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save plot: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = gRNAApp(root)
    root.mainloop()
