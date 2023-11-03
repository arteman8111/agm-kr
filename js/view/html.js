const render = (arr, label, index, title1 = `Заголовок`) => {
    const title = `<h2 style="margin: 25px 0 0 0; font-size: 24px; padding:10px 5px; ">${title1}</h3>`
    const app = document.querySelector('#app')
    const box = document.createElement('div')
    const table = document.createElement('table')
    table.className = "data"
    const header = document.createElement('tr')
    const td = (element) => `<td class="data__element">${element}</td>`
    header.insertAdjacentHTML('beforeend',td(''))
    for(let i = index; i < arr[0].length + index; i++){
        header.insertAdjacentHTML('beforeend', td(i))
    }
    for (let i = 0; i < arr.length; i++) {
        const tr = document.createElement('tr');
        tr.insertAdjacentHTML('beforeend',td(label[i]))
        arr[i].forEach(element => tr.insertAdjacentHTML('beforeend', td(math.format(element, {notation: 'exponential', precision: 5}))))
        table.insertAdjacentElement('beforeend', tr)
    }
    table.insertAdjacentElement('afterbegin', header)
    box.insertAdjacentHTML('beforeend', title)
    box.insertAdjacentElement('beforeend', table)
    app.insertAdjacentElement('beforeend', box)
}

export { render }