#ifndef PANGENOME_INDEX_CHUNKED_MEMORY_MAPPED_INDEX_HPP
#define PANGENOME_INDEX_CHUNKED_MEMORY_MAPPED_INDEX_HPP

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory>
#include <sstream>
#include "r-index.hpp"
#include "tag_arrays.hpp"

namespace panindexer {

/**
 * Chunked memory mapping for massive files
 * Maps files in manageable chunks to avoid virtual memory issues
 */
class ChunkedMemoryMapper {
private:
    int fd;
    size_t file_size;
    size_t chunk_size;
    std::vector<std::pair<void*, size_t>> mapped_chunks;
    bool is_open;
    
public:
    ChunkedMemoryMapper(size_t chunk_size_mb = 1024) 
        : fd(-1), file_size(0), chunk_size(chunk_size_mb * 1024 * 1024), is_open(false) {}
    
    ~ChunkedMemoryMapper() {
        close();
    }
    
    /**
     * Open file and prepare for chunked mapping
     * @param filename Path to file
     * @return true if successful
     */
    bool open(const std::string& filename) {
        if (is_open) {
            close();
        }
        
        fd = ::open(filename.c_str(), O_RDONLY);
        if (fd == -1) {
            std::cerr << "Cannot open file: " << filename << std::endl;
            return false;
        }
        
        struct stat st;
        if (fstat(fd, &st) == -1) {
            std::cerr << "Cannot get file stats for: " << filename << std::endl;
            ::close(fd);
            fd = -1;
            return false;
        }
        
        file_size = st.st_size;
        is_open = true;
        
        std::cerr << "Opened file for chunked mapping: " 
                  << (file_size / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
        
        return true;
    }
    
    /**
     * Map a specific chunk of the file
     * @param offset Starting offset in file
     * @param size Size of chunk to map
     * @return Pointer to mapped memory, or nullptr if failed
     */
    void* map_chunk(size_t offset, size_t size) {
        if (!is_open || fd == -1) {
            return nullptr;
        }
        
        // Ensure offset and size are within file bounds
        if (offset >= file_size) {
            return nullptr;
        }
        
        size_t actual_size = std::min(size, file_size - offset);
        
        // Map the chunk
        void* chunk = mmap(nullptr, actual_size, PROT_READ, MAP_PRIVATE, fd, offset);
        if (chunk == MAP_FAILED) {
            std::cerr << "Failed to map chunk at offset " << offset << std::endl;
            return nullptr;
        }
        
        mapped_chunks.emplace_back(chunk, actual_size);
        return chunk;
    }
    
    /**
     * Map the entire file in chunks
     * @return true if successful
     */
    bool map_all_chunks() {
        if (!is_open) return false;
        
        size_t offset = 0;
        while (offset < file_size) {
            size_t current_chunk_size = std::min(chunk_size, file_size - offset);
            
            void* chunk = map_chunk(offset, current_chunk_size);
            if (!chunk) {
                std::cerr << "Failed to map chunk at offset " << offset << std::endl;
                return false;
            }
            
            offset += current_chunk_size;
        }
        
        std::cerr << "Mapped " << mapped_chunks.size() << " chunks" << std::endl;
        return true;
    }
    
    /**
     * Get a specific chunk by index
     * @param chunk_index Index of chunk
     * @return Pair of (pointer, size) or (nullptr, 0) if invalid
     */
    std::pair<void*, size_t> get_chunk(size_t chunk_index) const {
        if (chunk_index >= mapped_chunks.size()) {
            return {nullptr, 0};
        }
        return mapped_chunks[chunk_index];
    }
    
    /**
     * Get total number of chunks
     */
    size_t get_chunk_count() const {
        return mapped_chunks.size();
    }
    
    /**
     * Close all mappings and file
     */
    void close() {
        // Unmap all chunks
        for (auto& [chunk, size] : mapped_chunks) {
            if (chunk != MAP_FAILED) {
                munmap(chunk, size);
            }
        }
        mapped_chunks.clear();
        
        if (fd != -1) {
            ::close(fd);
            fd = -1;
        }
        
        file_size = 0;
        is_open = false;
    }
    
    /**
     * Check if file is open
     */
    bool is_open_ready() const {
        return is_open;
    }
    
    /**
     * Get file size
     */
    size_t size() const {
        return file_size;
    }
    
    // Deleted copy constructor and assignment
    ChunkedMemoryMapper(const ChunkedMemoryMapper&) = delete;
    ChunkedMemoryMapper& operator=(const ChunkedMemoryMapper&) = delete;
};

/**
 * Chunked memory-mapped wrapper for FastLocate r-index
 * Maps large r-index files in chunks and deserializes them
 */
class ChunkedMemoryMappedRIndex {
private:
    ChunkedMemoryMapper mapper;
    FastLocate r_index;
    bool is_loaded;
    
public:
    ChunkedMemoryMappedRIndex(size_t chunk_size_mb = 1024) 
        : mapper(chunk_size_mb), is_loaded(false) {}
    
    /**
     * Open and load r-index file using chunked mapping
     * @param filename Path to .ri file
     * @return true if successful
     */
    bool open(const std::string& filename) {
        if (!mapper.open(filename)) {
            return false;
        }
        
        try {
            // For r-index files, we need to read the entire file
            // but we'll do it in chunks to avoid memory issues
            std::stringstream combined_stream;
            
            // Map all chunks and combine them
            if (!mapper.map_all_chunks()) {
                std::cerr << "Failed to map all chunks" << std::endl;
                return false;
            }
            
            // Combine all chunks into a single stream
            for (size_t i = 0; i < mapper.get_chunk_count(); ++i) {
                auto [chunk, size] = mapper.get_chunk(i);
                if (chunk) {
                    combined_stream.write(static_cast<const char*>(chunk), size);
                }
            }
            
            // Reset stream position and deserialize
            combined_stream.seekg(0);
            r_index.load(combined_stream);
            
            is_loaded = true;
            
            std::cerr << "Successfully loaded r-index using chunked memory mapping" << std::endl;
            return true;
            
        } catch (const std::exception& e) {
            std::cerr << "Failed to deserialize r-index: " << e.what() << std::endl;
            return false;
        }
    }
    
    /**
     * Close the mapper
     */
    void close() {
        mapper.close();
        is_loaded = false;
    }
    
    /**
     * Get reference to the loaded FastLocate object
     */
    FastLocate& get() {
        return r_index;
    }
    
    /**
     * Get const reference to the loaded FastLocate object
     */
    const FastLocate& get() const {
        return r_index;
    }
    
    /**
     * Check if index is loaded
     */
    bool is_ready() const {
        return is_loaded;
    }
    
    /**
     * Get file size
     */
    size_t size() const {
        return mapper.size();
    }
    
    ~ChunkedMemoryMappedRIndex() {
        close();
    }
};

/**
 * Chunked memory-mapped wrapper for TagArray
 * Maps large tag array files in chunks and deserializes them
 */
class ChunkedMemoryMappedTagArray {
private:
    ChunkedMemoryMapper mapper;
    TagArray tag_array;
    bool is_loaded;
    
public:
    ChunkedMemoryMappedTagArray(size_t chunk_size_mb = 1024) 
        : mapper(chunk_size_mb), is_loaded(false) {}
    
    /**
     * Open and load tag array file using chunked mapping
     * @param filename Path to .tags file
     * @return true if successful
     */
    bool open(const std::string& filename) {
        if (!mapper.open(filename)) {
            return false;
        }
        
        try {
            // Map all chunks
            if (!mapper.map_all_chunks()) {
                std::cerr << "Failed to map all chunks" << std::endl;
                return false;
            }
            
            // Combine chunks and deserialize
            std::stringstream combined_stream;
            
            for (size_t i = 0; i < mapper.get_chunk_count(); ++i) {
                auto [chunk, size] = mapper.get_chunk(i);
                if (chunk) {
                    combined_stream.write(static_cast<const char*>(chunk), size);
                }
            }
            
            combined_stream.seekg(0);
            tag_array.load_compressed_tags(combined_stream);
            
            is_loaded = true;
            
            std::cerr << "Successfully loaded tag array using chunked memory mapping" << std::endl;
            return true;
            
        } catch (const std::exception& e) {
            std::cerr << "Failed to deserialize tag array: " << e.what() << std::endl;
            return false;
        }
    }
    
    /**
     * Close the mapper
     */
    void close() {
        mapper.close();
        is_loaded = false;
    }
    
    /**
     * Get reference to the loaded TagArray object
     */
    TagArray& get() {
        return tag_array;
    }
    
    /**
     * Get const reference to the loaded TagArray object
     */
    const TagArray& get() const {
        return tag_array;
    }
    
    /**
     * Check if array is loaded
     */
    bool is_ready() const {
        return is_loaded;
    }
    
    /**
     * Get file size
     */
    size_t size() const {
        return mapper.size();
    }
    
    ~ChunkedMemoryMappedTagArray() {
        close();
    }
};

/**
 * Combined chunked memory-mapped index manager
 * Handles both r-index and tag array files efficiently using chunked mapping
 */
class ChunkedMemoryMappedIndexManager {
private:
    ChunkedMemoryMappedRIndex r_index;
    ChunkedMemoryMappedTagArray tag_array;
    size_t chunk_size_mb;
    
public:
    ChunkedMemoryMappedIndexManager(size_t chunk_size_mb = 1024) 
        : r_index(chunk_size_mb), tag_array(chunk_size_mb), chunk_size_mb(chunk_size_mb) {}
    
    /**
     * Open both r-index and tag array files using chunked mapping
     * @param r_index_file Path to .ri file
     * @param tag_array_file Path to .tags file
     * @return true if both files opened successfully
     */
    bool open(const std::string& r_index_file, const std::string& tag_array_file) {
        std::cerr << "Opening chunked memory-mapped indices (chunk size: " << chunk_size_mb << " MB)..." << std::endl;
        
        if (!r_index.open(r_index_file)) {
            std::cerr << "Failed to open r-index file" << std::endl;
            return false;
        }
        
        if (!tag_array.open(tag_array_file)) {
            std::cerr << "Failed to open tag array file" << std::endl;
            r_index.close();
            return false;
        }
        
        std::cerr << "Successfully opened both chunked memory-mapped indices" << std::endl;
        return true;
    }
    
    /**
     * Close both indices
     */
    void close() {
        r_index.close();
        tag_array.close();
    }
    
    /**
     * Get r-index reference
     */
    FastLocate& get_r_index() {
        return r_index.get();
    }
    
    /**
     * Get tag array reference
     */
    TagArray& get_tag_array() {
        return tag_array.get();
    }
    
    /**
     * Check if both indices are ready
     */
    bool is_ready() const {
        return r_index.is_ready() && tag_array.is_ready();
    }
    
    /**
     * Get total file size
     */
    size_t total_size() const {
        return r_index.size() + tag_array.size();
    }
    
    /**
     * Get chunk size in MB
     */
    size_t get_chunk_size_mb() const {
        return chunk_size_mb;
    }
    
    ~ChunkedMemoryMappedIndexManager() {
        close();
    }
};

} // namespace panindexer

#endif // PANGENOME_INDEX_CHUNKED_MEMORY_MAPPED_INDEX_HPP
