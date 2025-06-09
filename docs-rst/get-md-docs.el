;;; get-md-docs.el --- Extract Markdown documentation from source files.

;; Copyright (C) 2016 Mosè Giordano

;; get-md-docs is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 3, or (at your option)
;; any later version.

;; get-md-docs is distributed in the hope that it will be useful, but
;; WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;; General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with get-md-docs; see the file COPYING.  If not, write to the
;; Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
;; MA 02110-1301, USA.

;;; Code:

(defvar get-md-docs-files '("utils" "misc")
  "List of Markdown files to be created, without extension.")

(defvar get-md-docs--docs-regexp "\"\"\"\\(\\(?:.\\|\n\\)*\\)\"\"\""
  "Regexp matching Markdown docstrings in Julia source files.")

(defvar get-md-docs--file-name-regexp "^include(\"\\([^\"\n]+\\)\")$"
  "Regexp matching names of files included in the project.")

(defun get-md-docs-on-file (file-name)
  "Extract Markdown docstring in file FILE-NAME."
  (with-current-buffer
      (find-file-noselect file-name t)
    (goto-char (point-min))
    ;; Fetch the docstring using `get-md-docs--docs-regexp'.
    (when (re-search-forward get-md-docs--docs-regexp nil t)
      (replace-regexp-in-string
       ;; Replace escaped backslashes (needed in Markdown) with single
       ;; backslashes.
       "\\\\\\\\" "\\\\"
       (replace-regexp-in-string
	;; Replace "\$" with "$" for math modes.
	"\\\\\\$" "\$"
	(replace-regexp-in-string
	 ;; Demote sections for the online docs, from "###" to "####".
	 "^### \\([^#\n]*\\) ###$"
	 "##### \\1 ####"
	 (match-string-no-properties 1)))))))

(defun get-md-docs-search-file (file-name)
  "Search in FILE-NAME files from which extract docstrings.

The real extraction is then done by `get-md-docs-on-file'."
  (let ((backup-inhibited backup-inhibited)
	file docs)
    ;; Visit the Julia file that includes the source of each function.
    (with-current-buffer
	(find-file-noselect (concat "../src/" file-name ".jl") t)
      (goto-char (point-min))
      (save-excursion
	(save-match-data
	  ;; Find the names of Julia files that load each function and then
	  ;; operate on them.
	  (while (re-search-forward get-md-docs--file-name-regexp nil t)
	    (setq
	     file (match-string-no-properties 1)
	     docs
	     ;; Append the new docstring to the previous ones in order to
	     ;; collect them all.
	     (concat
	      docs
	      ;; Add a header with the name of the function.
	      "\n### " (file-name-base file) " ###\n"
	      ;; Extract the docstring from `file'.
	      (get-md-docs-on-file file)
	      ;; Add a horizontal rule to clearly separate sections.
	      "\n----\n"))))))
    ;; After collecting the docstrings, write the result to the Markdown file.
    (with-current-buffer
	(find-file-noselect (concat file-name ".md") t)
      (setq backup-inhibited t)
      (erase-buffer)
      (insert docs)
      ;; Remove the last horizontal rule.
      (if (re-search-backward "----" nil t)
	  (replace-match ""))
      (delete-trailing-whitespace)
      (save-buffer))))

(defun get-md-docs ()
  "Extract all Mardown docstring and save to file."
  (dolist (file get-md-docs-files)
    (get-md-docs-search-file file)))

(defvar fix-rst-docs--function-regexp "^~+\n\\(\n+::[ \t\n\r]+\\)"
  "Regexp matching the function definition.")

(defun fix-rst-doc-file (file-name)
  "Fix formatting of rst FILE-NAME."
  (let ((backup-inhibited t)
	file docs)
    (with-current-buffer
	(find-file-noselect (concat file-name ".rst") t)
      (goto-char (point-min))
      (save-excursion
	(save-match-data
	  ;; Turn the beginning of docstring of each function info a `function'
	  ;; block.
	  (while (re-search-forward
		  (concat "\\(?:"
			  ;; Group 1.
			  fix-rst-docs--function-regexp
			  "\\|"
			  ;; Group 2.
			  "\\(^\.\. code::\\)"
			  "\\)")
		  nil t)
	    (cond
	     ((match-beginning 1)
	      ;; Turn the beginning of docstring of each function info a
	      ;; `function' block.
	      (replace-match "\n.. function:: " nil nil nil 1))
	     ;; Replace ".. code::" with ".. code-block::"
	     ((match-beginning 2)
	      (replace-match ".. code-block::"))))))
      (delete-trailing-whitespace)
      (save-buffer))))

(defun fix-rst-docs ()
  "Run formatting improvements to every rst files."
  (dolist (file get-md-docs-files)
    (fix-rst-doc-file file)))

;;; get-md-docs.el ends here
