#!/usr/bin/env python3
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
import subprocess
import threading

class MatrixEncryptGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Matrix Encryption Tool")
        self.root.geometry("700x600")
        self.root.resizable(True, True)
        
        self.setup_styles()
        
        self.binary_path = self.find_binary()
        self.setup_ui()
    
    def setup_styles(self):
        style = ttk.Style()
        style.configure('Title.TLabel', font=('Arial', 16, 'bold'))
        style.configure('Status.TLabel', font=('Arial', 10))
        style.configure('Execute.TButton', font=('Arial', 12, 'bold'))
    
    def find_binary(self):
        possible_names = ['./matrixencrypt', 'matrixencrypt', './matrixencrypt.exe']
        for name in possible_names:
            if os.path.exists(name):
                return name
        return None
    
    def setup_ui(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Title
        title_label = ttk.Label(main_frame, text="Matrix Encryption Tool", style='Title.TLabel')
        title_label.pack(pady=(0, 20))
        
        # Binary status
        self.setup_binary_status(main_frame)
        
        # Mode selection
        self.setup_mode_selection(main_frame)
        
        # File selection
        self.setup_file_selection(main_frame)
        
        # Execute button
        self.setup_execute_section(main_frame)
        
        # Log area
        self.setup_log_area(main_frame)
    
    def setup_binary_status(self, parent):
        status_frame = ttk.Frame(parent)
        status_frame.pack(fill=tk.X, pady=(0, 10))
        
        binary_text = "Binary status: FOUND" if self.binary_path else "Binary status: NOT FOUND"
        color = "green" if self.binary_path else "red"
        
        self.binary_label = ttk.Label(status_frame, text=binary_text, foreground=color)
        self.binary_label.pack(side=tk.LEFT)
        
        if not self.binary_path:
            compile_btn = ttk.Button(status_frame, text="Compile C++ Code", command=self.compile_cpp)
            compile_btn.pack(side=tk.RIGHT)
    
    def setup_mode_selection(self, parent):
        mode_frame = ttk.LabelFrame(parent, text="Operation Mode", padding="10")
        mode_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.mode_var = tk.StringVar(value="encrypt")
        
        ttk.Radiobutton(mode_frame, text="🔒 Encrypt File", 
                       variable=self.mode_var, value="encrypt").pack(side=tk.LEFT, padx=20)
        ttk.Radiobutton(mode_frame, text="🔓 Decrypt File", 
                       variable=self.mode_var, value="decrypt").pack(side=tk.LEFT, padx=20)
    
    def setup_file_selection(self, parent):
        files_frame = ttk.Frame(parent)
        files_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Key file
        self.key_file = tk.StringVar()
        self.create_file_row(files_frame, "Key File:", self.key_file, 0)
        
        # Input file  
        self.input_file = tk.StringVar()
        self.create_file_row(files_frame, "Input File:", self.input_file, 1)
        
        # Output file
        self.output_file = tk.StringVar()
        self.create_file_row(files_frame, "Output File:", self.output_file, 2, is_output=True)
    
    def create_file_row(self, parent, label, variable, row, is_output=False):
        row_frame = ttk.Frame(parent)
        row_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(row_frame, text=label, width=12).pack(side=tk.LEFT)
        
        entry = ttk.Entry(row_frame, textvariable=variable, width=50)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        
        browse_text = "Save As..." if is_output else "Browse"
        browse_cmd = lambda: self.browse_output_file(variable) if is_output else self.browse_input_file(variable)
        ttk.Button(row_frame, text=browse_text, command=browse_cmd).pack(side=tk.RIGHT)
    
    def browse_input_file(self, variable):
        filename = filedialog.askopenfilename(
            title="Select Input File",
            filetypes=[("All files", "*.*"), ("Text files", "*.txt")]
        )
        if filename:
            variable.set(filename)
    
    def browse_output_file(self, variable):
        filename = filedialog.asksaveasfilename(
            title="Save Output File As",
            defaultextension=".txt",
            filetypes=[("All files", "*.*"), ("Text files", "*.txt")]
        )
        if filename:
            variable.set(filename)
    
    def setup_execute_section(self, parent):
        execute_frame = ttk.Frame(parent)
        execute_frame.pack(fill=tk.X, pady=10)
        
        self.execute_btn = ttk.Button(execute_frame, text="🚀 Execute Operation", 
                                     command=self.execute_operation, style='Execute.TButton')
        self.execute_btn.pack(pady=5)
        
        # Progress
        self.progress = ttk.Progressbar(execute_frame, mode='indeterminate')
        
        # Status
        self.status_var = tk.StringVar(value="Ready")
        status_label = ttk.Label(execute_frame, textvariable=self.status_var, style='Status.TLabel')
        status_label.pack(pady=5)
    
    def setup_log_area(self, parent):
        log_frame = ttk.LabelFrame(parent, text="Execution Log", padding="10")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=(10, 0))
        
        self.log_text = scrolledtext.ScrolledText(log_frame, height=15, width=80)
        self.log_text.pack(fill=tk.BOTH, expand=True)
        
        # Add copy functionality
        self.log_text.bind("<Control-c>", self.copy_log)
    
    def copy_log(self, event=None):
        try:
            selected = self.log_text.get(tk.SEL_FIRST, tk.SEL_LAST)
            self.root.clipboard_clear()
            self.root.clipboard_append(selected)
        except tk.TclError:
            # Nothing selected
            pass
        return "break"
    
    def log_message(self, message):
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()
    
    def compile_cpp(self):
        if not os.path.exists('main.cpp'):
            messagebox.showerror("Error", "main.cpp not found in current directory!")
            return
        
        self.log_message("compiling C++ code...")
        self.execute_btn.config(state='disabled')
        
        def compile_thread():
            try:
                result = subprocess.run(
                    ['g++', '-o', 'matrixencrypt', 'main.cpp', '-std=c++17', '-O2'],
                    capture_output=True, text=True
                )
                
                self.root.after(0, self.compile_finished, result)
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Error", f"Compilation failed: {str(e)}"))
        
        threading.Thread(target=compile_thread, daemon=True).start()
    
    def compile_finished(self, result):
        self.execute_btn.config(state='normal')
        
        if result.returncode == 0:
            self.log_message("Compilation successful!")
            self.binary_path = './matrixencrypt'
            self.binary_label.config(text="Binary status: FOUND", foreground="green")
            messagebox.showinfo("Success", "C++ code compiled successfully!")
        else:
            self.log_message("Compilation failed!")
            self.log_message(result.stderr)
            messagebox.showerror("Compilation Error", result.stderr)
    
    def validate_inputs(self):
        if not self.binary_path:
            messagebox.showerror("Error", "C++ binary not found!\nPlease compile the code first.")
            return False
        
        required = {
            "Key file": self.key_file.get(),
            "Input file": self.input_file.get(), 
            "Output file": self.output_file.get()
        }
        
        missing = [name for name, value in required.items() if not value.strip()]
        if missing:
            messagebox.showerror("Error", f"Please specify:\n" + "\n".join(f"• {name}" for name in missing))
            return False
        
        # Check if input files exist
        for file_type, file_path in [("Key file", self.key_file.get()), 
                                   ("Input file", self.input_file.get())]:
            if not os.path.exists(file_path):
                messagebox.showerror("Error", f"{file_type} not found:\n{file_path}")
                return False
        
        return True
    
    def execute_operation(self):
        if not self.validate_inputs():
            return
        
        # Disable button and show progress
        self.execute_btn.config(state='disabled')
        self.progress.pack(fill=tk.X, pady=5)
        self.progress.start()
        self.status_var.set("Executing...")
        
        mode = self.mode_var.get()
        
        def operation_thread():
            try:
                mode_flag = "-encrypt" if mode == "encrypt" else "-decrypt"
                
                command = [
                    self.binary_path,
                    mode_flag,
                    self.key_file.get(),
                    "-o", 
                    self.output_file.get(),
                    self.input_file.get()
                ]
                
                self.root.after(0, lambda: self.log_message(f"💻 Command: {' '.join(command)}"))
                
                result = subprocess.run(command, capture_output=True, text=True)
                
                self.root.after(0, self.operation_finished, result)
                
            except Exception as e:
                self.root.after(0, lambda: self.operation_finished(None, str(e)))
        
        threading.Thread(target=operation_thread, daemon=True).start()
    
    def operation_finished(self, result, exception=None):
        # Stop progress and re-enable button
        self.progress.stop()
        self.progress.pack_forget()
        self.execute_btn.config(state='normal')
        
        if exception:
            self.status_var.set("Error")
            self.log_message(f"Exception: {exception}")
            messagebox.showerror("Error", f"Operation failed:\n{exception}")
            return
        
        if result.returncode == 0:
            self.status_var.set("Success!")
            self.log_message("Operation completed successfully!")
            self.log_message(f"Output saved to: {self.output_file.get()}")
            messagebox.showinfo("Success", 
                              f"Operation completed!\n\nOutput file: {self.output_file.get()}")
        else:
            self.status_var.set("Error")
            error_msg = result.stderr if result.stderr else "Unknown error occurred"
            self.log_message(f"Error: {error_msg}")
            messagebox.showerror("Error", f"Operation failed:\n{error_msg}")

def main():
    # Check if tkinter is available
    try:
        root = tk.Tk()
    except ImportError:
        print("ERROR: tkinter not available!")
        print("Please install tkinter:")
        return 1
    
    app = MatrixEncryptGUI(root)
    root.mainloop()
    return 0

if __name__ == "__main__":
    exit(main())