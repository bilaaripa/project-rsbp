// download.js

function downloadImage(imagePath) {
  // Buat elemen <a> secara dinamis
  var link = document.createElement("a");
  link.href = imagePath;
  link.download = "dendrogram.png";

  // Simulasikan klik pada tautan
  document.body.appendChild(link);
  link.click();

  // Hapus elemen setelah klik
  document.body.removeChild(link);
}
