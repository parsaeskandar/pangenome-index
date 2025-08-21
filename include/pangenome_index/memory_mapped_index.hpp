#ifndef PANGENOME_INDEX_MEMORY_MAPPED_INDEX_HPP
#define PANGENOME_INDEX_MEMORY_MAPPED_INDEX_HPP

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include "r-index.hpp"
#include "tag_arrays.hpp"

namespace panindexer {

/**
 * Memory-mapped wrapper for FastLocate r-index
 * Uses memory mapping for efficient file access while properly deserializing SDSL structures
 */
class MemoryMappedRIndex {
private:
    int fd;
    void* mapped_data;
    size_t file_size;
    FastLocate r_index;
    bool is_open;
    
public:
    MemoryMappedRIndex() : fd(-1), mapped_data(MAP_FAILED), file_size(0), is_open(false) {}
    
    ~MemoryMappedRIndex() {
        close();
    }
    
    /**
     * Open and memory-map an r-index file, then deserialize it
     * @param filename Path to the .ri file
     * @return true if successful, false otherwise
     */
    bool open(const std::string& filename) {
        if (is_open) {
            close();
        }
        
        // Open the file
        fd = ::open(filename.c_str(), O_RDONLY);
        if (fd == -1) {
            std::cerr << "Cannot open r-index file: " << filename << std::endl;
            return false;
        }
        
        // Get file size
        struct stat st;
        if (fstat(fd, &st) == -1) {
            std::cerr << "Cannot get file stats for: " << filename << std::endl;
            ::close(fd);
            fd = -1;
            return false;
        }
        
        file_size = st.st_size;
        
        // Memory map the file
        mapped_data = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapped_data == MAP_FAILED) {
            std::cerr << "Memory mapping failed for: " << filename << std::endl;
            ::close(fd);
            fd = -1;
            return false;
        }
        
        // Create a string stream from the mapped memory and deserialize
        try {
            std::string mapped_string(static_cast<const char*>(mapped_data), file_size);
            std::istringstream stream(mapped_string);
            
            // Use the existing load method to properly deserialize
            r_index.load(stream);
            
        } catch (const std::exception& e) {
            std::cerr << "Failed to deserialize r-index: " << e.what() << std::endl;
            close();
            return false;
        }
        
        is_open = true;
        
        std::cerr << "Successfully memory-mapped and deserialized r-index: " 
                  << (file_size / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
        
        return true;
    }
    
    /**
     * Close the memory mapping and file
     */
    void close() {
        if (mapped_data != MAP_FAILED) {
            munmap(mapped_data, file_size);
            mapped_data = MAP_FAILED;
        }
        if (fd != -1) {
            ::close(fd);
            fd = -1;
        }
        file_size = 0;
        is_open = false;
    }
    
    /**
     * Get reference to the deserialized FastLocate object
     * @return Reference to FastLocate
     */
    FastLocate& get() {
        return r_index;
    }
    
    /**
     * Get const reference to the deserialized FastLocate object
     * @return Const reference to FastLocate
     */
    const FastLocate& get() const {
        return r_index;
    }
    
    /**
     * Check if the index is open and ready
     */
    bool is_open_ready() const {
        return is_open;
    }
    
    /**
     * Get the file size in bytes
     */
    size_t size() const {
        return file_size;
    }
    
    // Deleted copy constructor and assignment to prevent accidental copying
    MemoryMappedRIndex(const MemoryMappedRIndex&) = delete;
    MemoryMappedRIndex& operator=(const MemoryMappedRIndex&) = delete;
};

/**
 * Memory-mapped wrapper for TagArray
 * Uses memory mapping for efficient file access while properly deserializing compressed data
 */
class MemoryMappedTagArray {
private:
    int fd;
    void* mapped_data;
    size_t file_size;
    TagArray tag_array;
    bool is_open;
    
public:
    MemoryMappedTagArray() : fd(-1), mapped_data(MAP_FAILED), file_size(0), is_open(false) {}
    
    ~MemoryMappedTagArray() {
        close();
    }
    
    /**
     * Open and memory-map a tag array file, then deserialize it
     * @param filename Path to the .tags file
     * @return true if successful, false otherwise
     */
    bool open(const std::string& filename) {
        if (is_open) {
            close();
        }
        
        // Open the file
        fd = ::open(filename.c_str(), O_RDONLY);
        if (fd == -1) {
            std::cerr << "Cannot open tag array file: " << filename << std::endl;
            return false;
        }
        
        // Get file size
        struct stat st;
        if (fstat(fd, &st) == -1) {
            std::cerr << "Cannot get file stats for: " << filename << std::endl;
            ::close(fd);
            fd = -1;
            return false;
        }
        
        file_size = st.st_size;
        
        // Memory map the file
        mapped_data = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapped_data == MAP_FAILED) {
            std::cerr << "Memory mapping failed for: " << filename << std::endl;
            ::close(fd);
            fd = -1;
            return false;
        }
        
        // Create a string stream from the mapped memory and deserialize
        try {
            std::string mapped_string(static_cast<const char*>(mapped_data), file_size);
            std::istringstream stream(mapped_string);
            
            // Use the existing load_compressed_tags method to properly deserialize
            tag_array.load_compressed_tags(stream);
            
        } catch (const std::exception& e) {
            std::cerr << "Failed to deserialize tag array: " << e.what() << std::endl;
            close();
            return false;
        }
        
        is_open = true;
        
        std::cerr << "Successfully memory-mapped and deserialized tag array: " 
                  << (file_size / (1024.0 * 1024.0 * 1024.0)) << " GB" << std::endl;
        
        return true;
    }
    
    /**
     * Close the memory mapping and file
     */
    void close() {
        if (mapped_data != MAP_FAILED) {
            munmap(mapped_data, file_size);
            mapped_data = MAP_FAILED;
        }
        if (fd != -1) {
            ::close(fd);
            fd = -1;
        }
        file_size = 0;
        is_open = false;
    }
    
    /**
     * Get reference to the deserialized TagArray object
     * @return Reference to TagArray
     */
    TagArray& get() {
        return tag_array;
    }
    
    /**
     * Get const reference to the deserialized TagArray object
     * @return Const reference to TagArray
     */
    const TagArray& get() const {
        return tag_array;
    }
    
    /**
     * Check if the array is open and ready
     */
    bool is_open_ready() const {
        return is_open;
    }
    
    /**
     * Get the file size in bytes
     */
    size_t size() const {
        return file_size;
    }
    
    // Deleted copy constructor and assignment to prevent accidental copying
    MemoryMappedTagArray(const MemoryMappedTagArray&) = delete;
    MemoryMappedTagArray& operator=(const MemoryMappedTagArray&) = delete;
};

/**
 * Combined memory-mapped index manager
 * Handles both r-index and tag array files efficiently
 */
class MemoryMappedIndexManager {
private:
    MemoryMappedRIndex r_index;
    MemoryMappedTagArray tag_array;
    
public:
    MemoryMappedIndexManager() = default;
    
    /**
     * Open both r-index and tag array files
     * @param r_index_file Path to .ri file
     * @param tag_array_file Path to .tags file
     * @return true if both files opened successfully
     */
    bool open(const std::string& r_index_file, const std::string& tag_array_file) {
        std::cerr << "Opening memory-mapped indices..." << std::endl;
        
        if (!r_index.open(r_index_file)) {
            std::cerr << "Failed to open r-index file" << std::endl;
            return false;
        }
        
        if (!tag_array.open(tag_array_file)) {
            std::cerr << "Failed to open tag array file" << std::endl;
            r_index.close();
            return false;
        }
        
        std::cerr << "Successfully opened both memory-mapped indices" << std::endl;
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
        return r_index.is_open_ready() && tag_array.is_open_ready();
    }
    
    /**
     * Get total memory-mapped size
     */
    size_t total_size() const {
        return r_index.size() + tag_array.size();
    }
    
    ~MemoryMappedIndexManager() {
        close();
    }
};

} // namespace panindexer

#endif // PANGENOME_INDEX_MEMORY_MAPPED_INDEX_HPP
