# Makefile

.PHONY: deploy-docs

# just a simple file to copy the outputs from the sphinx build into the docs dir.. this is all for github pages
DOCS_DIR := docs
SOURCE_BUILD_DIR := docs_maker/build/html

deploy-docs:
	@echo "--- Starting documentation deployment ---"

	# Remove all contents from the target 'docs' directory
	@echo "Cleaning contents of $(DOCS_DIR)/..."
	rm -rf $(DOCS_DIR)
	mkdir $(DOCS_DIR)


	# Copy all files from the source build directory to 'docs'
	@echo "Copying built HTML from $(SOURCE_BUILD_DIR)/ to $(DOCS_DIR)/..."
	cp -r $(SOURCE_BUILD_DIR)/. $(DOCS_DIR)/

	@echo "--- Documentation deployment complete ---"